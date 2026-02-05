# BN Thickness Estimation Tools (Interactive Matplotlib)

This folder contains small Python scripts for extracting color values along user-selected lines/ROIs on an optical micrograph and (optionally) mapping the **normalized Red-channel contrast** to **hBN thickness** after calibration.

---

## Requirements

- Python 3.9+ recommended
- Packages:
  - `numpy`
  - `matplotlib`
  - `pillow`
  - `scipy` (only needed for spline-based thickness calibration/inversion)

Install (pip):

```bash
pip install numpy matplotlib pillow scipy
```

> **Linux note:** `tkinter` is needed for the file picker.
> - Ubuntu/Debian: `sudo apt-get install python3-tk`

---

## Scripts Overview

### 1) Line Profile GUI (R-only normalized)

This script lets you:
- Select a **background** region (for `R_bg`)
- Draw **multiple lines**
- Plot **normalized R** along each line:
  - If background set: `R_norm = (R - R_bg) / R_bg`
  - Otherwise: `R_norm = R / 255`
- Export CSV for each line and a combined CSV

#### Run
```bash
python thickness_check_multi.py
```

When it starts, a file picker opens. Choose an image (PNG/JPG/TIF/...).

#### Controls (typical)
- **Click 2 points**: create a new line (active line)
- **Drag endpoints**: adjust the active line
- **Enter**: finalize the current line and start a new one
- **Set BG (click)**: click a clean substrate area to define background color
- **Save CSV / Save ALL**: export profiles
- **Close window**: finish

> Button names may vary slightly depending on your version of the script, but the behavior above is the intended workflow.

#### Output
- Per-line CSVs: e.g. `line_profile_01.csv`, `line_profile_02.csv`, ...
- Combined CSV: e.g. `line_profiles_ALL.csv`

Typical columns include:
- `dist_px` (distance along line in pixels)
- `R` (0–255)
- `R_norm` (normalized)
- `x`, `y` sample coordinates  
(Your script may include additional columns if enabled.)

---

### 2) Calibration Tool (ROI → Thickness Model)

This script helps build a **calibration curve** between **thickness (nm)** and **R_norm** using ROIs:
1. Drag a rectangle on clean substrate → defines `R_bg`
2. Drag rectangles on flakes with known thickness labels → enter thickness when prompted
3. Fits a spline `R_norm(t)` and saves model + points

#### Run
```bash
python calibrate_bn_thickness_from_Rnorm.py
```

Select your calibration image when prompted.

#### Output
- `calibration_points.csv` — measured `(t_nm, R_norm)` points
- `calibration_model.npz` — background and spline parameters

---

## Using Calibration in Your GUI

Once you have `calibration_points.csv` (or `.npz`), your GUI can:
1. Compute `R_norm` for a measurement
2. Estimate thickness by inverting the fitted curve

Recommended approach:
- Use a **spline / interpolant** (robust)
- Inversion may return **multiple candidates** (see limitations)

---

## Limitations / Important Notes

1. **Non-unique thickness mapping**
   - Optical contrast vs thickness on SiO₂/Si is generally **oscillatory** (thin-film interference). This works for BN on 285 nm SiO2 and only for range 10-30 ish nm BN.
   - The same `R_norm` may correspond to **multiple thicknesses**.
   - Best practice: constrain the thickness range, or use **multiple channels** (R/G/B or luminance) to disambiguate.

2. **Calibration is microscope-specific**
   - White balance, illumination spectrum, camera response, exposure, polarizers, objective NA, etc. change the mapping.
   - Always calibrate with the **same microscope settings** used for unknown samples.

3. **Background selection matters**
   - Pick background on clean substrate near the flake, avoiding dust/edges/gradients.
   - The scripts typically use a **median** over a patch to reduce noise, but strong gradients still affect results.

4. **Image processing assumptions**
   - Sampling along lines may use **nearest** pixel (default) or **bilinear** interpolation (if enabled).
   - JPEG compression and aggressive denoising/sharpening can distort channel values.

5. **Units**
   - Line distance is reported in **pixels** by default.
   - Convert to µm using a known scale bar (future enhancement: click scale bar endpoints).

6. **GUI backend**
   - File picker uses `tkinter`. If the picker does not open, install `python3-tk` (Linux) or run in a desktop environment.

---

## Quick Troubleshooting

- **No file picker / crashes on import tkinter**
  - Install Tk:
    - Ubuntu/Debian: `sudo apt-get install python3-tk`

- **Dragging endpoints doesn't work**
  - Make sure you click near the endpoint marker.
  - Some Matplotlib backends behave differently; try running from a standard desktop Python (not inside certain IDE plot panes).

- **Weird normalization**
  - Ensure you set background (`Set BG`) before trusting `(R - R_bg)/R_bg`.

---

## License / Notes

These scripts are intended for research/analysis workflows and can be modified freely for your lab use.
# Flake_thickness_Estimation
