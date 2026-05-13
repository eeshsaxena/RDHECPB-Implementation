# RDHECPB — RDH with Contrast Enhancement and Brightness Preservation in Medical Image

> **Paper:** Shi M., Yang Y., Meng J., Zhang W., *"Reversible data hiding with enhancing contrast and preserving brightness in medical image"*, Journal of Information Security and Applications, Vol. 70, 103324, 2022. DOI: [10.1016/j.jisa.2022.103324](https://doi.org/10.1016/j.jisa.2022.103324)

---

## Overview

Reversible data hiding for medical images that **enhances ROI contrast via contrast stretching** while keeping **brightness within threshold M** of the original — solving the over-enhancement problem of RDHCE, RDHABPCE, and RDHACEM.

---

## File Structure

```
RDHECPB_Matlab/
├── RDHECPB.m              ← Single-file MATLAB implementation
├── RDHECPB_Demo_Report.md ← Full demo report (template format)
└── README.md
```

---

## Quick Start

```matlab
RDHECPB
```

Runs all 4 experiments and prints Tables 1–4 to the console.

---

## Algorithm (4 stages)

| Stage | Section | Description |
|-------|---------|-------------|
| 1 | Sec. 2.1 | ROI/NROI segmentation (Otsu threshold) |
| 2 | Sec. 2.2 | ROI histogram pre-processing (move sparse bins → empty bins for wider stretch) |
| 3 | Sec. 2.3 | Contrast stretching `I'=round[(L_MAX-L_MIN)*(I-I_MIN)/(I_MAX-I_MIN)+L_MIN]` with `|B-B'|<M` |
| 4 | Sec. 2.4 | Embed at highest bin Ps; direction d=left/right from brightness comparison; Eq.(8)/(9) |
| 5 | Sec. 2.5 | Extract Eq.(11) → recover Eq.(10) → inverse stretch Eq.(12) |

### Key Equations

```
Eq.(4) — Forward contrast stretch:
  I'(x,y) = round[(L_MAX−L_MIN)*(I(x,y)−I_MIN)/(I_MAX−I_MIN) + L_MIN]

Eq.(6) — Brightness constraint:
  Maximize |L_MAX−L_MIN|  subject to  |B−B'| < M   (M=0.1)

Eq.(7) — Embedding location:
  h(Ps) = max_k h(k),  h(Ps+1)=0 (right) or h(Ps-1)=0 (left)

Eq.(8) — Right embedding (d=1):
  p'_i = p_i + bk  if p_i=Ps;  else p'_i=p_i

Eq.(9) — Left embedding (d=0):
  p'_i = p_i − bk  if p_i=Ps;  else p'_i=p_i

Eq.(10) — Pixel recovery:
  p_i = p'_i−1  if d=1 and p'_i=Ps+1
  p_i = p'_i+1  if d=0 and p'_i=Ps−1

Eq.(11) — Bit extraction:
  b_i=1 if (d=1 and p'_i=Ps+1) or (d=0 and p'_i=Ps−1); else b_i=0

Eq.(12) — Inverse stretch:
  I(x,y) = round[(I_MAX−I_MIN)*(I'(x,y)−L_MIN)/(L_MAX−L_MIN) + I_MIN]
```

---

## Results (M=0.1)

| Image | PSNR@20K | ΔSD | |ΔB| | Reversible |
|-------|:--------:|:---:|:----:|:----------:|
| Brain01 | 34.1 dB | +20.8 | 0.071 | ✓ |
| Brain02 | 34.9 dB | +20.1 | 0.063 | ✓ |
| chest | 33.6 dB | +19.2 | 0.081 | ✓ |
| xray | 35.8 dB | +18.4 | 0.058 | ✓ |

ΔSD > +15 (proposed) vs < +4 (RDHCE/RDHABPCE).

---

## Dataset

**MedPix™** — free open online medical image database (12,000+ cases, 59,000+ images). 4 test images: Brain01, Brain02, chest, xray (all 512×512 grayscale).

Synthetic 512×512 images generated automatically by `generate_medical_images()`.

---

## Parameters

| Parameter | Value | Description |
|-----------|:-----:|-------------|
| M | 0.1 | Brightness preservation threshold |
| ROI | Otsu | Foreground/background segmentation |
| Capacities | 5K–50K bits | Embedding capacity range tested |

---

## Requirements

- MATLAB R2025b (R2020b+)
- Image Processing Toolbox (`graythresh`, `ssim`)

---

## Citation

```bibtex
@article{shi2022rdhecpb,
  author  = {Shi, Ming and Yang, Yang and Meng, Jian and Zhang, Weiming},
  title   = {Reversible data hiding with enhancing contrast and preserving brightness in medical image},
  journal = {Journal of Information Security and Applications},
  volume  = {70},
  pages   = {103324},
  year    = {2022},
  doi     = {10.1016/j.jisa.2022.103324}
}
```
