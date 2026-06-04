# Speckle-based Digital Image Correlation (DIC)
### GRISHMA Undergraduate Summer Research Program — IIT Kharagpur
**Project:** Development of Microstructure-based Digital Image Correlation (Micro-DIC) in Small-Scale Mechanical Testing
**Faculty:** Prof. Shibayan Roy | Materials Science Center
**Project Ref:** GRISHMA-FAC-2026-0140

---

## Table of Contents
1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Python Version](#python-version)
4. [MATLAB Version](#matlab-version)
5. [Algorithm Description](#algorithm-description)
6. [Parameters Guide](#parameters-guide)
7. [Output Description](#output-description)
8. [Results Summary](#results-summary)
9. [Known Limitations](#known-limitations)
10. [Development Log](#development-log)

---

## Overview

This repository contains two implementations of a **Speckle-based Digital Image Correlation (DIC)** pipeline developed for small-scale mechanical testing. DIC is a non-contact optical technique that tracks a random speckle pattern applied to a specimen surface to measure full-field displacement and strain during deformation.

### What DIC does
```
Reference image (unloaded)  →  DIC Algorithm  →  Displacement field (U, V)
Deformed image (loaded)                        →  Strain field (εxx, εyy, εxy)
                                               →  Von Mises strain
                                               →  Poisson ratio
```

### Key Features
- **Sub-pixel accuracy** via parabolic interpolation (MATLAB) / ICGN refinement (Python)
- **Grip-referenced displacement** — left grip = zero, measures true elongation
- **LSQ strain calculation** — least-squares polynomial fit over a window (more accurate than finite differences)
- **Outlier rejection** — MAD-based detection and replacement
- **Automatic masking** — low-correlation points masked and shown in red
- **Publication-quality plots** — white background, colorbars, titles

---

## Repository Structure

```
DIC_Project/
│
├── Python/
│   └── dic_ncorr.py          ← Main Python DIC script (Ncorr-style)
│
├── MATLAB/
│   └── dic_main.m            ← Main MATLAB DIC script
│
├── README.md                 ← This file
└── DIC_SOP.docx             ← Standard Operating Procedure
```

---

## Python Version

### Requirements

```bash
pip install numpy opencv-python matplotlib scipy
```

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥ 1.21 | Array operations |
| opencv-python | ≥ 4.5 | Image loading, NCC template matching |
| matplotlib | ≥ 3.4 | Plotting and figure saving |
| scipy | ≥ 1.7 | Spline interpolation, Gaussian filtering |

### Setup

1. Place all images in a folder (e.g. `D:\Images\`)
2. Open `dic_ncorr.py` in any editor
3. Edit the USER SETTINGS section at the top
4. Run: `python dic_ncorr.py`

### Configuration (`dic_ncorr.py`)

```python
# ── USER SETTINGS ──────────────────────────────────────────

IMAGE_FOLDER   = r"D:\Micro DIC Project\Speckle Pattern DIC\Images"
IMAGE_PATTERN  = "*.JPG"
RESULTS_FOLDER = r"D:\Micro DIC Project\Speckle Pattern DIC\dic_results_final"

# ROI — gauge section only (x1, y1, x2, y2) in pixels
ROI = (973, 853, 1596, 982)

# DIC parameters
SUBSET_SIZE    = 25      # correlation window size (pixels, odd number)
STEP_SIZE      = 4       # grid spacing (pixels)
MAX_SHIFT      = 25      # coarse search range (pixels)
ICGN_ITER      = 10      # ICGN refinement iterations
ICGN_TOL       = 1e-4    # ICGN convergence tolerance

# Post-processing
STRAIN_WINDOW  = 7       # LSQ strain window (grid points, odd)
STRAIN_SMOOTH  = 1.2     # Gaussian smoothing sigma
OUTLIER_NSIGMA = 2.5     # outlier rejection threshold
MIN_CORR       = 0.6     # minimum correlation threshold
REF_GRIP_COLS  = 3       # columns used as fixed reference

# Physical scale (optional)
GAUGE_LENGTH_MM = None   # set to gauge length in mm if known
```

### How to Find Your ROI

Open `img01.JPG` in any image viewer. Find the pixel coordinates of the four corners of the **speckle-covered gauge section only** (not the grips). Use those as `(x1, y1, x2, y2)`.

In Python/OpenCV you can also do:
```python
import cv2
img = cv2.imread('img01.JPG')
cv2.imshow('Select ROI', img)
roi = cv2.selectROI('Select ROI', img)  # drag to select, press Enter
print(roi)  # prints (x, y, w, h) → convert to (x, x+w, y, y+h)
```

### Running

```bash
python dic_ncorr.py
```

Expected terminal output:
```
======================================================================
  Ncorr-style DIC — GRISHMA Project, IIT Kharagpur
  ICGN refinement | LSQ strain | Grip-referenced
======================================================================

Found      : 78 images
Reference  : img01.JPG
ROI        : (973, 853, 1596, 982)  →  623 × 129 px

Frame 001/077  img02.JPG   ZNCC=0.9868  max_U=+0.24px  exx=+0.022%
Frame 002/077  img03.JPG   ZNCC=0.9902  max_U=+1.05px  exx=+0.163%  nu=0.48
...
```

### Python Method — ICGN (Inverse Compositional Gauss-Newton)

The Python version uses a two-stage approach:

**Stage 1 — Coarse NCC search:**
- Template matching with `cv2.matchTemplate`
- Parabolic sub-pixel interpolation at correlation peak
- Provides initial displacement estimate (u₀, v₀)

**Stage 2 — ICGN refinement:**
- Bicubic spline interpolation of both images (5th order)
- Iteratively minimizes the ZNCC residual
- Converges to sub-pixel accuracy (typically < 0.05 px error)
- Falls back to coarse result if refinement diverges

---

## MATLAB Version

### Requirements

- MATLAB R2022a or later (R2026a recommended)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

Verify with: `ver`

### Setup (MATLAB Online)

1. Go to [matlab.mathworks.com](https://matlab.mathworks.com)
2. Log in with your institutional account
3. Upload all 78 images to MATLAB Drive
4. Upload `dic_main.m` to MATLAB Drive
5. Open `dic_main.m` in the editor
6. Edit USER SETTINGS
7. Press F5 or click Run

### Configuration (`dic_main.m`)

```matlab
%% ── USER SETTINGS ───────────────────────────────────────────

IMAGE_FOLDER   = '/MATLAB Drive';
IMAGE_PATTERN  = 'img*.JPG';
RESULTS_FOLDER = '/MATLAB Drive/dic_results';

% ROI — gauge section only (x1, y1, x2, y2)
ROI = [973, 853, 1596, 982];

% DIC parameters
SUBSET_SIZE  = 29;     % correlation window (pixels, odd)
STEP_SIZE    = 5;      % grid spacing (pixels)
MAX_SHIFT    = 40;     % coarse search range (pixels)

% Post-processing
STRAIN_WINDOW  = 5;    % LSQ strain window (grid points)
STRAIN_SMOOTH  = 1.2;  % Gaussian smoothing sigma
OUTLIER_NSIGMA = 2.5;  % outlier rejection threshold
MIN_CORR       = 0.70; % minimum correlation threshold
REF_GRIP_COLS  = 3;    % columns used as fixed reference

% Physical scale (optional)
GAUGE_LENGTH_MM = [];  % set to gauge length in mm if known
```

### How to Find Your ROI in MATLAB

```matlab
img = imread('img01.JPG');
imshow(img);
roi = round(getPosition(imrect()));  % drag rectangle, double-click to confirm
% roi = [x, y, width, height]
% Convert to: ROI = [x, y, x+width, y+height]
```

Or use Paint: hover cursor over corners of gauge section,
read coordinates at bottom of screen.

### Running

Open `dic_main.m` → Press F5

Expected terminal output:
```
=======================================================
  DIC Pipeline — GRISHMA Project, IIT Kharagpur
=======================================================
Found      : 78 images
Reference  : img01.JPG
ROI        : [973 853 1596 982] → 623 × 129 px
Subset     : 29px | Step: 5px | Max shift: 40px
Grid size  : 110 × 11 = 1210 points

Frame 001/077  img02.JPG   ZNCC=0.9878  maxU=+0.23px  exx=+0.020%
Frame 002/077  img03.JPG   ZNCC=0.9910  maxU=+1.07px  exx=+0.166%  nu=0.26
...
```

### MATLAB Method — NCC + Parabolic Sub-pixel

The MATLAB version uses:

**Stage 1 — NCC search:**
- `normxcorr2` for normalized cross-correlation
- Finds integer-pixel peak in correlation map
- Parabolic interpolation for sub-pixel correction

**Stage 2 — LSQ strain:**
- Fits a linear polynomial over a local window of displacement values
- More accurate than simple finite differences
- Handles NaN points from masked regions automatically

---

## Algorithm Description

### Step 1 — Image Preprocessing
```
Load image → Convert to grayscale → Normalize (zero mean, unit variance)
```
Normalization makes the correlation insensitive to brightness changes between frames.

### Step 2 — Grid Generation
```
Half = SUBSET_SIZE // 2
Valid region = [Half+MAX_SHIFT : W-Half-MAX_SHIFT] × [Half+MAX_SHIFT : H-Half-MAX_SHIFT]
Grid points spaced STEP_SIZE pixels apart
```

### Step 3 — Correlation
For each grid point (cx, cy):
```
1. Extract reference subset: ref_img[cy-half:cy+half, cx-half:cx+half]
2. Search window in deformed image: ±MAX_SHIFT around (cx, cy)
3. NCC: find peak location → integer displacement (dx_int, dy_int)
4. Sub-pixel: fit parabola → fractional displacement (sub_dx, sub_dy)
5. Total displacement: dx = dx_int + sub_dx
```

### Step 4 — Post-processing
```
Outlier rejection (MAD-based) → Grip reference subtraction → Gaussian smoothing
```

### Step 5 — Strain Calculation (LSQ Window)
For each grid point (i, j):
```
Collect displacement values in window of size STRAIN_WINDOW × STRAIN_WINDOW
Fit: u(x,y) = a + b·x + c·y  (linear polynomial)
εxx = du/dx = b / STEP_SIZE
εyy = dv/dy = d / STEP_SIZE  (from v polynomial)
εxy = 0.5 × (du/dy + dv/dx)
```

### Step 6 — Output
```
Displacement maps (U, V) → Strain maps (εxx, εyy, εxy, Von Mises) → Plots saved as PNG
```

---

## Parameters Guide

| Parameter | Effect | Typical Range |
|-----------|--------|---------------|
| `SUBSET_SIZE` | Larger = more stable but lower spatial resolution | 21–61 px |
| `STEP_SIZE` | Smaller = denser grid, more detail, slower | 3–10 px |
| `MAX_SHIFT` | Must be ≥ maximum displacement between consecutive frames | 15–60 px |
| `STRAIN_WINDOW` | Larger = smoother strain, loses localization | 3–9 grid pts |
| `STRAIN_SMOOTH` | Post-smoothing of strain field | 0.5–2.5 |
| `MIN_CORR` | Threshold below which points are masked | 0.60–0.90 |
| `OUTLIER_NSIGMA` | Lower = more aggressive outlier removal | 2.0–4.0 |
| `REF_GRIP_COLS` | Number of reference columns for zero-point | 2–5 |

### Choosing SUBSET_SIZE
```
Rule of thumb:
- Speckle diameter ≈ 3–5 px → SUBSET_SIZE = 21–29
- Speckle diameter ≈ 5–10 px → SUBSET_SIZE = 29–41
- Must satisfy: SUBSET_SIZE + MAX_SHIFT < ROI_height / 2
```

### Choosing MAX_SHIFT
```
Check your displacement data:
- If max elongation is 50 px over 78 frames
- Between consecutive frames: ~50/77 ≈ 0.65 px (very small)
- But camera shake adds ±5–10 px
- So MAX_SHIFT = 15–25 is usually enough
- If frames show exactly MAX_SHIFT×2 displacement → increase MAX_SHIFT
```

---

## Output Description

### Frame Plots (`frame_001.png` to `frame_077.png`)

Each frame plot contains 8 panels:

| Panel | Content |
|-------|---------|
| ROI + Grid | Reference image with grid overlay. Green = valid, Red = masked |
| ZNCC Correlation | Spatial map of correlation quality. Green > 0.95 = excellent |
| U displacement | Horizontal displacement relative to left grip (pixels) |
| V displacement | Vertical displacement relative to left grip (pixels) |
| Von Mises strain | Combined equivalent strain |
| εxx (axial) | Strain in loading direction (positive = tension) |
| εyy (transverse) | Strain perpendicular to loading (negative = Poisson contraction) |
| εxy (shear) | Shear strain (should be near zero for pure tension) |

### Summary Plot (`summary.png`)

| Panel | Content |
|-------|---------|
| Axial strain evolution | Mean and 95th percentile εxx vs frame |
| Transverse strain evolution | Mean εyy vs frame |
| Gauge elongation | Max U displacement vs frame |
| Poisson ratio evolution | Apparent ν vs frame with median line |

### Data Files

**Python:** `data_001.npz` to `data_077.npz`
```python
import numpy as np
data = np.load('data_001.npz')
ux   = data['ux']        # horizontal displacement field
uy   = data['uy']        # vertical displacement field
exx  = data['exx']       # axial strain field
eyy  = data['eyy']       # transverse strain field
exy  = data['exy']       # shear strain field
corr = data['corr']      # correlation coefficient field
mask = data['bad_mask']  # boolean mask (True = bad point)
```

**MATLAB:** `data_001.mat` to `data_077.mat`
```matlab
load('data_001.mat')
% Variables: grid_x, grid_y, ux, uy, exx, eyy, exy, corr, bad_mask
```

---

## Results Summary

Results obtained from 78-image tensile test sequence:

| Parameter | Value |
|-----------|-------|
| Specimen | Unknown metal (awaiting confirmation from Prof. Roy) |
| Gauge ROI | 623 × 129 px |
| Grid points | 1210 (110 × 11) |
| Mean ZNCC (elastic region) | 0.987–0.991 |
| Poisson ratio ν | 0.391 (median, frames 8–26) |
| Peak axial strain | ~7.4–7.8% |
| Max gauge elongation | ~50 px |
| Fracture frame | ~38–39 |
| Necking onset | ~Frame 27 |

### Loading Stages

| Frames | Stage | ZNCC | Notes |
|--------|-------|------|-------|
| 1–17 | Elastic loading | > 0.985 | Linear strain increase |
| 17–27 | Plastic deformation | 0.976–0.985 | ν begins rising |
| 27–38 | Necking | 0.957–0.976 | Strain localizing at right end |
| 38–53 | Post-necking | 0.920–0.957 | Large deformation |
| 54–70 | Post-fracture | 0.78–0.83 | Unreliable — specimen broken |
| 71–77 | Unloaded | 0.82–0.85 | Elastic springback |

---

## Known Limitations

1. **MAX_SHIFT capping** — At very large deformations (frame 37–38), displacement reaches exactly MAX_SHIFT limit. Results may be slightly underestimated. Increase MAX_SHIFT to 60 if needed.

2. **Thin gauge section** — ROI is only 129px tall (11 grid rows). This limits the resolution of transverse strain measurement and makes ε_yy noisier than ε_xx.

3. **No physical scale** — Gauge length in mm is not yet confirmed. All displacement values are in pixels. To convert: `mm = px × (gauge_length_mm / gauge_length_px)`.

4. **Reference-based DIC only** — All frames are correlated against frame 1. At very high strains (>15%), correlation degrades. Incremental DIC (frame-to-frame) would be more accurate at large deformations.

5. **Out-of-plane motion** — V displacement shows a slight gradient from top to bottom, suggesting minor out-of-plane rotation of the specimen. This introduces a small error in ε_yy.

6. **No machine data validation** — Load-displacement data from the testing machine has not yet been obtained for cross-validation.

---

## Development Log

| Version | Changes |
|---------|---------|
| v1 | Basic NCC with integer-pixel displacement; noisy strain maps |
| v2 | Added parabolic sub-pixel interpolation; reduced noise |
| v3 | Corrected ROI to gauge section only; added rigid body correction |
| v4 | Grip-referenced displacement; Poisson ratio recovered correctly |
| v5 (Python) | ICGN refinement added; bicubic spline interpolation |
| v5 (MATLAB) | NCC + LSQ strain; white background plots; MATLAB Online compatible |
| Current | Both versions stable; publication-quality output |

---

## Contact

**Student:** GRISHMA-STU-2026-0363
**Faculty:** Prof. Shibayan Roy (shibayan@matsc.iitkgp.ac.in)
**Institution:** IIT Kharagpur, Materials Science Center
**Program:** GRISHMA 2026 | KRITI Initiative

---

*This code was developed as part of the GRISHMA Undergraduate Summer Research Program 2026 at IIT Kharagpur. The speckle-based DIC pipeline is the first phase of a larger project to develop Microstructure-based DIC (Micro-DIC) for small-scale mechanical testing.*
