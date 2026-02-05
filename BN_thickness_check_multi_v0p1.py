import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.widgets import Button
from tkinter import Tk
from tkinter.filedialog import askopenfilename
# Xuetao Ma


# ---------------------- helpers ----------------------
def estimate_bn_thickness_from_norm_rgb(step_norm: float) -> float:
    """
    Estimate BN thickness (nm) from a *background-normalized* R value:

        norm_rgb = (R - R_bg) / R_bg

    Using your background (from the screenshot):
        BG RGB = [160, 127, 140]  ->  R_bg = 160

    Calibration data (raw R measured on flake, assumed linear):
        R_raw           = [50, 82, 140, 102]
        AFM_thickness   = [17, 22, 27, 23]   (nm)

    We convert R_raw -> norm_rgb using R_bg=160, then fit:
        thickness_nm = m * norm_rgb + b
    """
    import numpy as np
    R_bg = 160.0  # BG RGB = [160,127,140] -> use R channel

    R_diff = np.array([50.0, 82.0, 140.0, 102.0], dtype=float)
    t_nm  = np.array([17.0, 22.0, 27.0, 23.0], dtype=float)

    x = R_diff / R_bg
    X = np.column_stack([x, np.ones_like(x)])
    m, b = np.linalg.lstsq(X, t_nm, rcond=None)[0]
    return float(m * float(step_norm) + b)

def pick_image_file():
    """Open a native file dialog and return selected image path (or None if canceled)."""
    root = Tk()
    root.withdraw()           # hide main window
    root.attributes("-topmost", True)
    path = askopenfilename(
        title="Select an image",
        filetypes=[
            ("Image files", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.webp"),
            ("All files", "*.*"),
        ],
    )
    root.destroy()
    return path if path else None

def luminance_bt601(rgb):
    """
    rgb: (N,3) array, dtype uint8 or float
    returns: (N,) float luminance in same scale as rgb (0-255 if rgb is 0-255)
    """
    rgb = rgb.astype(float)
    R, G, B = rgb[:, 0], rgb[:, 1], rgb[:, 2]
    return 0.299 * R + 0.587 * G + 0.114 * B
    
def sample_rgb_along_line(img_rgb, p0, p1, n_samples=800, method="nearest"):
    """
    Returns:
      dist_px: (N,) float
      rgb: (N,3) uint8/float
      xs, ys: (N,) float sample coordinates
    """
    x0, y0 = p0
    x1, y1 = p1
    xs = np.linspace(x0, x1, n_samples)
    ys = np.linspace(y0, y1, n_samples)

    H, W = img_rgb.shape[:2]

    if method == "nearest":
        xi = np.clip(np.rint(xs).astype(int), 0, W - 1)
        yi = np.clip(np.rint(ys).astype(int), 0, H - 1)
        rgb = img_rgb[yi, xi, :].astype(np.float32)
    elif method == "bilinear":
        # bilinear interpolation for smoother profiles
        x0i = np.floor(xs).astype(int)
        y0i = np.floor(ys).astype(int)
        x1i = np.clip(x0i + 1, 0, W - 1)
        y1i = np.clip(y0i + 1, 0, H - 1)
        x0i = np.clip(x0i, 0, W - 1)
        y0i = np.clip(y0i, 0, H - 1)

        dx = (xs - x0i)[:, None]
        dy = (ys - y0i)[:, None]

        Ia = img_rgb[y0i, x0i].astype(np.float32)
        Ib = img_rgb[y0i, x1i].astype(np.float32)
        Ic = img_rgb[y1i, x0i].astype(np.float32)
        Id = img_rgb[y1i, x1i].astype(np.float32)

        rgb = (Ia * (1 - dx) * (1 - dy) +
               Ib * dx       * (1 - dy) +
               Ic * (1 - dx) * dy +
               Id * dx       * dy)
    else:
        raise ValueError("method must be 'nearest' or 'bilinear'")

    dist_px = np.linspace(0, np.hypot(x1 - x0, y1 - y0), n_samples)
    return dist_px, rgb, xs, ys

def plateau_values_from_step(y, guard=30, smooth_win=9):
    """
    Find largest step in y and return (left_plateau, right_plateau, step_index).

    guard: number of points to ignore near ends & near the step
    smooth_win: moving average window for robust derivative (odd recommended)
    """
    y = np.asarray(y, float)
    n = y.size
    if n < 2 * guard + 10:
        return float(np.median(y)), float(np.median(y)), None

    # Smooth a bit for robust step detection
    w = int(max(3, smooth_win))
    if w % 2 == 0:
        w += 1
    kernel = np.ones(w) / w
    ys = np.convolve(y, kernel, mode="same")

    # Find index of maximum absolute gradient away from edges
    dy = np.diff(ys)
    lo = guard
    hi = max(lo + 1, n - guard - 1)
    idx_local = np.argmax(np.abs(dy[lo:hi]))
    i0 = lo + idx_local  # step between i0 and i0+1

    # Take medians on both sides, away from the step by 'guard'
    left_end = max(guard, i0 - guard)
    right_start = min(n - guard, i0 + 1 + guard)

    left = np.median(y[guard:left_end]) if left_end > guard else np.median(y[:i0+1])
    right = np.median(y[right_start:n-guard]) if (n-guard) > right_start else np.median(y[i0+1:])

    return float(left), float(right), int(i0)

def robust_background_rgb(img_rgb, x, y, half_window=10):
    """
    Estimate background RGB from a (2*half_window+1)^2 square around (x,y)
    using the median (robust to dust / flakes).
    """
    H, W = img_rgb.shape[:2]
    xi = int(np.clip(round(x), 0, W - 1))
    yi = int(np.clip(round(y), 0, H - 1))

    x0 = max(0, xi - half_window)
    x1 = min(W, xi + half_window + 1)
    y0 = max(0, yi - half_window)
    y1 = min(H, yi + half_window + 1)

    patch = img_rgb[y0:y1, x0:x1, :].reshape(-1, 3).astype(np.float32)
    bg = np.median(patch, axis=0)  # (3,)
    return bg


# ---------------------- config ----------------------
# If you run from the same folder as the image, keep just filename.
# Or give an absolute path.
img_path = None

N_SAMPLES = 800
METHOD = "nearest"   # "nearest" or "bilinear"
BG_HALF_WINDOW = 12  # background median over (2*win+1)^2 around click

# ---------------------- load image ----------------------
if img_path is None:
    img_path = pick_image_file()
    if img_path is None:
        raise SystemExit("No image selected. Exiting.")

#img_rgb = np.array(Image.open(IMG_PATH).convert("RGB"))
img_rgb = np.array(Image.open(img_path).convert("RGB"))

# ---------------------- state ----------------------
state = {
    "mode": "pick_line_p0",   # pick_line_p0 -> pick_line_p1 -> active_line (draggable) -> pick_line_p0 ...
    "bg_mode": False,         # if True, next click sets background
    "bg_rgb": None,           # (3,) float
    "plot_mode": "delta",     # "delta" or "raw"
    "active": None,           # dict for current line
    "lines": [],              # list of finalized line dicts
}

# line dict structure:
# {
#   "p0": (x,y), "p1": (x,y),
#   "line_artist": Line2D,
#   "p0_artist": PathCollection, "p1_artist": PathCollection,
#   "dist": (N,), "rgb": (N,3), "xs":(N,), "ys":(N,)
# }

# ---------------------- figure ----------------------
fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0])

ax_img = fig.add_subplot(gs[0, 0])
ax_prof = fig.add_subplot(gs[0, 1])

ax_img.imshow(img_rgb)
ax_img.axis("off")

# Profile lines (single set = for active line)
r_line, = ax_prof.plot([], [], label="R (normalized)")

ax_prof.set_xlabel("Distance (pixels)")
ax_prof.set_ylabel("Value")
ax_prof.set_title("RGB profile (active line)")
ax_prof.legend(loc="best")

# Buttons
btn_save_ax = fig.add_axes([0.80, 0.02, 0.17, 0.06])
btn_save = Button(btn_save_ax, "Save ALL CSV")

btn_bg_ax = fig.add_axes([0.60, 0.02, 0.18, 0.06])
btn_bg = Button(btn_bg_ax, "Set BG (click)")

btn_mode_ax = fig.add_axes([0.40, 0.02, 0.18, 0.06])
btn_mode = Button(btn_mode_ax, "Plot: Delta")

status_text = fig.text(0.01, 0.01, "", fontsize=10)

def set_status(msg):
    status_text.set_text(msg)
    fig.canvas.draw_idle()

def update_title():
    if state["active"] is None:
        ax_prof.set_title("R profile (no active line)")
        return

    line = state["active"]
    if "rgb" not in line:
        ax_prof.set_title("R profile (active line)")
        return

    rgb = line["rgb"].astype(np.float32)
    R = rgb[:, 0]  # raw R (0..255)

    # If BG is set, use it; else infer BG from the line using the step
    if state.get("bg_rgb", None) is not None:
        R_bg = float(state["bg_rgb"][0])
        # Use step to infer flake plateau only (more robust than overall median)
        left, right, i0 = plateau_values_from_step(R)
        # pick flake plateau as the one farthest from BG
        R_flake = left if abs(left - R_bg) > abs(right - R_bg) else right
        R_norm = (R_flake - R_bg) / max(R_bg, 1e-9)
        mode_str = "step vs BG"
    else:
        # No BG: infer both plateaus and define BG as the higher plateau (often substrate brighter in R)
        left, right, i0 = plateau_values_from_step(R)
        R_bg = max(left, right)
        R_flake = min(left, right)
        R_norm = (R_flake - R_bg) / max(R_bg, 1e-9)
        mode_str = "step (auto BG)"
        
    step_norm = abs(left - right) / max(R_bg, 1e-9)
    t_est = estimate_bn_thickness_from_norm_rgb(step_norm)

    ax_prof.set_title(
        f"R profile | {mode_str} | R_norm(step)={R_norm:+.4f} → thickness≈{t_est:.1f} nm"
    )

def compute_and_store(line_dict):
    p0, p1 = line_dict["p0"], line_dict["p1"]
    dist, rgb, xs, ys = sample_rgb_along_line(img_rgb, p0, p1, n_samples=N_SAMPLES, method=METHOD)
    line_dict["dist"] = dist
    line_dict["rgb"] = rgb
    line_dict["xs"] = xs
    line_dict["ys"] = ys

def refresh_profile_from_active():
    if state["active"] is None:
        r_line.set_data([], [])
        ax_prof.set_xlim(0, 1)
        ax_prof.set_ylim(0, 1)
        ax_prof.set_ylabel("R (normalized)")
        fig.canvas.draw_idle()
        return

    line = state["active"]
    compute_and_store(line)

    dist = line["dist"]
    rgb = line["rgb"].astype(np.float32)
    R = rgb[:, 0]

    if state["bg_rgb"] is not None:
        R_bg = float(state["bg_rgb"][0])
        R_norm = (R - R_bg) / max(R_bg, 1e-9)
        ax_prof.set_ylabel("(R - R_bg) / R_bg")
    else:
        R_norm = R / 255.0
        ax_prof.set_ylabel("R / 255")

    r_line.set_data(dist, R_norm)
    ax_prof.set_xlim(dist.min(), max(dist.max(), 1))

    ymin, ymax = float(np.min(R_norm)), float(np.max(R_norm))
    pad = 0.05 * (ymax - ymin + 1e-9)
    ax_prof.set_ylim(ymin - pad, ymax + pad)

    update_title()
    fig.canvas.draw_idle()

def nearest_endpoint(x, y, p0, p1, thresh=12):
    d0 = np.hypot(x - p0[0], y - p0[1])
    d1 = np.hypot(x - p1[0], y - p1[1])
    if d0 <= thresh and d0 <= d1:
        return "p0"
    if d1 <= thresh:
        return "p1"
    return None

def finalize_active_line():
    """Freeze current active line (no longer draggable), keep on image, store it."""
    if state["active"] is None:
        return
    state["lines"].append(state["active"])
    state["active"] = None
    state["mode"] = "pick_line_p0"
    refresh_profile_from_active()
    set_status(f"Finalized a line. Total lines: {len(state['lines'])}. Now click to start a new line.")

def start_new_active(p0):
    # create fresh artists for the new active line
    line_artist, = ax_img.plot([], [], linewidth=2)
    p0_artist = ax_img.scatter([], [], s=60)
    p1_artist = ax_img.scatter([], [], s=60)

    state["active"] = {
        "p0": p0,
        "p1": p0,
        "line_artist": line_artist,
        "p0_artist": p0_artist,
        "p1_artist": p1_artist,
    }
    state["mode"] = "pick_line_p1"

def update_active_artists():
    if state["active"] is None:
        return
    p0 = state["active"]["p0"]
    p1 = state["active"]["p1"]
    state["active"]["line_artist"].set_data([p0[0], p1[0]], [p0[1], p1[1]])
    state["active"]["p0_artist"].set_offsets([p0])
    state["active"]["p1_artist"].set_offsets([p1])

def on_click(event):
    if event.inaxes != ax_img or event.xdata is None or event.ydata is None:
        return
    x, y = float(event.xdata), float(event.ydata)

    # background setting mode
    if state["bg_mode"]:
        state["bg_rgb"] = robust_background_rgb(img_rgb, x, y, half_window=BG_HALF_WINDOW)
        state["bg_mode"] = False
        update_title()
        refresh_profile_from_active()
        set_status("Background set. Now draw lines (2 clicks per line).")
        return

    # left click: line creation / selecting endpoint drag
    if state["mode"] == "pick_line_p0":
        start_new_active((x, y))
        update_active_artists()
        update_title()
        set_status("Picked first point. Click second point to set the line end.")
        fig.canvas.draw_idle()
        return

    if state["mode"] == "pick_line_p1":
        state["active"]["p1"] = (x, y)
        state["mode"] = "active_line"
        update_active_artists()
        refresh_profile_from_active()
        set_status("Line created. Drag endpoints to adjust. Click 'Save ALL CSV' anytime. "
                   "To start another line: press Enter (finalize) or click 'Save' then Enter.")
        return

    if state["mode"] == "active_line":
        # initiate dragging if clicked near an endpoint
        p0 = state["active"]["p0"]
        p1 = state["active"]["p1"]
        which = nearest_endpoint(x, y, p0, p1, thresh=12)
        state["dragging"] = which
        return

def on_release(event):
    state["dragging"] = None

def on_move(event):
    if state.get("dragging") is None:
        return
    if state["active"] is None:
        return
    if event.inaxes != ax_img or event.xdata is None or event.ydata is None:
        return
    x, y = float(event.xdata), float(event.ydata)

    if state["dragging"] == "p0":
        state["active"]["p0"] = (x, y)
    elif state["dragging"] == "p1":
        state["active"]["p1"] = (x, y)

    update_active_artists()
    refresh_profile_from_active()

def on_key(event):
    # Enter finalizes active line and allows starting a new one
    if event.key in ("enter", "return"):
        if state["mode"] == "active_line":
            finalize_active_line()
            update_title()
    # Escape cancels current active line
    if event.key == "escape":
        if state["active"] is not None:
            # remove its artists from the axes
            state["active"]["line_artist"].remove()
            state["active"]["p0_artist"].remove()
            state["active"]["p1_artist"].remove()
            state["active"] = None
            state["mode"] = "pick_line_p0"
            refresh_profile_from_active()
            update_title()
            set_status("Canceled active line. Now click to start a new line.")

def on_save(event):
    """
    Saves ALL finalized lines + (optional) current active line into CSV files.
    """
    out_dir = os.path.abspath(".")
    base = "rgb_profiles"

    to_save = list(state["lines"])
    if state["active"] is not None and state["mode"] in ("pick_line_p1", "active_line"):
        # include active line as a preview (not added to finalized list)
        to_save = to_save + [state["active"]]

    if len(to_save) == 0:
        set_status("No lines to save yet.")
        return

    combined_rows = []
    for i, ln in enumerate(to_save, start=1):
        # ensure data exists
        if "dist" not in ln:
            compute_and_store(ln)

        dist = ln["dist"]
        rgb = ln["rgb"].astype(np.float32)
        xs = ln["xs"]
        ys = ln["ys"]

        if state["bg_rgb"] is None:
            bg = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        else:
            bg = state["bg_rgb"].astype(np.float32)

        delta = rgb - bg[None, :]

        # per-line file
        fname = f"{base}_line{i:02d}.csv"
        out = np.column_stack([dist, rgb[:,0], rgb[:,1], rgb[:,2],
                               delta[:,0], delta[:,1], delta[:,2],
                               xs, ys])
        header = "dist_px,R,G,B,dR,dG,dB,x,y"
        np.savetxt(fname, out, delimiter=",", header=header, comments="")

        # combined rows with line index
        combined_rows.append(np.column_stack([np.full_like(dist, i), out]))

    combined = np.vstack(combined_rows)
    combined_name = f"{base}_ALL.csv"
    np.savetxt(combined_name, combined, delimiter=",",
               header="line_id,dist_px,R,G,B,dR,dG,dB,x,y", comments="")

    set_status(f"Saved {len(to_save)} line(s): {base}_lineXX.csv + {combined_name} in {out_dir}")

def on_set_bg(event):
    state["bg_mode"] = True
    set_status("Background selection: click a BACKGROUND region (no flakes). Using median in a square patch.")

def on_toggle_plotmode(event):
    if state["plot_mode"] == "delta":
        state["plot_mode"] = "raw"
        btn_mode.label.set_text("Plot: Raw")
    else:
        state["plot_mode"] = "delta"
        btn_mode.label.set_text("Plot: Delta")
    refresh_profile_from_active()

btn_save.on_clicked(on_save)
btn_bg.on_clicked(on_set_bg)
btn_mode.on_clicked(on_toggle_plotmode)

fig.canvas.mpl_connect("button_press_event", on_click)
fig.canvas.mpl_connect("motion_notify_event", on_move)
fig.canvas.mpl_connect("button_release_event", on_release)
fig.canvas.mpl_connect("key_press_event", on_key)

update_title()
set_status("1) (Optional) Click 'Set BG' then click background.\n"
           "2) Click 2 points to draw a line. Drag endpoints to adjust.\n"
           "3) Press Enter to finalize the line and draw another. Close window to finish.")
plt.tight_layout()
plt.show()
