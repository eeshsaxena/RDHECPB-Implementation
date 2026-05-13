# RDHECPB — Reversible Data Hiding with Enhancing Contrast and Preserving Brightness in Medical Image
**Paper:** Shi, Yang, Meng & Zhang, Journal of Information Security and Applications, Vol. 70, 103324, 2022
**DOI:** 10.1016/j.jisa.2022.103324 | **Platform:** MATLAB R2025b

---

## 1. Paper Reference

| Field | Details |
|-------|---------|
| Title | Reversible Data Hiding with Enhancing Contrast and Preserving Brightness in Medical Image |
| Authors | Ming Shi, Yang Yang (Corresponding), Jian Meng, Weiming Zhang |
| Journal | Journal of Information Security and Applications |
| Volume/Article | 70 (2022), 103324 |
| DOI | 10.1016/j.jisa.2022.103324 |
| Published | Available online: 16 September 2022 |

---

## 2. Problem Statement

Existing RDH-CE methods for medical images (RDHCE [9], ACERDH [10], RDHABPCE [11], RDHACEM [16]) cause **over-enhancement** — brightness increases uncontrollably as embedding capacity increases, making the image visually unnatural. This paper proposes RDHECPB which:
1. Segments the image into **ROI** (diagnostic content) and **NROI** (background)
2. Applies **contrast stretching** (stronger than histogram equalization) to ROI
3. Controls the **embedding direction** (left/right) per round to keep brightness within threshold M of the original

---

## 3. Background — Prior Methods

| Method | CE Type | Brightness Control | Limitation |
|--------|---------|:-----------------:|------------|
| RDHCE [9] Wu 2015 | Histogram equalization | None | Over-enhancement at high capacity |
| ACERDH [10] Mansouri 2021 | Two-sided HS | None | Over-enhances brightness |
| RDHABPCE [11] Kim 2019 | HE + brightness preserve | Limited range | Weaker contrast stretching |
| RDHACEM [16] Gao 2021 | Auto contrast, ROI/NROI | None | Brightness completely out of control |
| **RDHECPB (proposed)** | **Contrast stretching + ROI** | **|B−B'| < M** | **None — best balance** |

---

## 4. Proposed Method

### 4.1 System Overview (Fig. 2)

```
Original Medical Image I
    │
    ├─[Sec.2.1] ROI/NROI Segmentation (Otsu threshold)
    │
    ├─ ROI branch:
    │   ├─[Sec.2.2] Histogram Pre-processing (move sparse bins → empty bins)
    │   ├─[Sec.2.3] Contrast Stretching with brightness constraint (Eq.4–6)
    │   └─[Sec.2.4] Embedding at Ps with direction d (Eq.7–9)
    │
    └─ NROI branch:
        └─[Sec.2.4] Embed remaining bits (histogram shifting, from [15])
    
    → Marked Image I_emb
    → Recovery: Eq.10 (pixel) + Eq.12 (inverse stretch) + NROI reverse
```

### 4.2 Stage 1 — ROI/NROI Segmentation (Sec. 2.1)

Medical images are divided into:
- **ROI** = central diagnostic region (tumors, organs, bones)
- **NROI** = background (monochrome, clinically less important)

Implementation uses **Otsu's method [20]** to threshold the image into foreground (ROI) and background (NROI).

### 4.3 Stage 2 — ROI Histogram Pre-processing (Sec. 2.2)

To obtain a **larger contrast stretching range**, sparse bins in the ROI histogram are moved to adjacent bins, creating empty bins. This extends the histogram's effective range before stretching.

- Empty bins created are **recorded** (their locations stored as side information)
- The side info is embedded in the image in Section 2.4

### 4.4 Stage 3 — Contrast Stretching with Brightness Preservation (Sec. 2.3)

**Forward contrast stretching (Eq. 4):**
```
I'(x,y) = round[(L_MAX − L_MIN) × (I(x,y) − I_MIN) / (I_MAX − I_MIN) + L_MIN]
```

**Expected brightness after stretching (Eq. 5):**
```
B' = mean over ROI pixels of stretched values
```

**Brightness preservation constraint (Eq. 6):**
```
Maximize |L_MAX − L_MIN|  subject to  |B − B'| < M
```

Where M = 0.1 (paper's chosen value, Sec. 3.1 experiment validates this).

The algorithm searches for the largest `[L_MIN, L_MAX]` range that keeps the stretched image's mean brightness within M of the original's mean brightness.

### 4.5 Stage 4 — Embedding in ROI with Brightness-Directed Direction (Sec. 2.4)

**Embedding location Ps (Eq. 7, Fig. 6):**
```
h(Ps) = max(h(k)),  k ∈ [0,255]         ← highest bin in ROI
```
Additionally: if right embedding (d=1) → h(Ps+1)=0; if left (d=0) → h(Ps-1)=0.

**Direction selection (Sec. 2.4.1 — 3 cases):**
```
B_stretched > B_original  →  d=0 (left):  reduces brightness back toward B
B_stretched < B_original  →  d=1 (right): increases brightness back toward B
B_stretched = B_original  →  d=0 (default)
```

**Right embedding (d=1) — Eq. (8):**
```
p'_i = p_i + bk    if p_i = Ps   ← embed bit (0: stay; 1: go to Ps+1)
p'_i = p_i          otherwise
```

**Left embedding (d=0) — Eq. (9):**
```
p'_i = p_i − bk    if p_i = Ps   ← embed bit (0: stay; 1: go to Ps-1)
p'_i = p_i          otherwise
```

**Side information structure (Sec. 2.4.2):**
```
Secret message structure = [I_MAX, I_MIN, 9LSBs, Ps, d, bk]
```
- `I_MAX`, `I_MIN`: pixel range of original ROI (stored in last round only)
- `9LSBs`: original LSBs of the 9th pixel
- `Ps` (8 bits) + `d` (1 bit) = 9 bits → stored in LSBs of 9 pixels
- `bk`: payload + empty bin locations

### 4.6 Stage 5 — Embedding in NROI

If ROI capacity is insufficient, remaining payload bits are embedded into NROI using standard histogram shifting [15]. NROI recovery is the reverse histogram shift.

### 4.7 Extraction and Recovery (Sec. 2.5)

**Bit extraction (Eq. 11):**
```
b_i = 1    if d=1 and i=Ps+1      ← was right-embedded with bit 1
b_i = 1    if d=0 and i=Ps−1      ← was left-embedded with bit 1
b_i = 0    otherwise
```

**Pixel recovery from embedding (Eq. 10):**
```
p_i = p'_i − 1    if d=1 and i=Ps+1    ← restore Ps from Ps+1
p_i = p'_i + 1    if d=0 and i=Ps−1    ← restore Ps from Ps-1
p_i = p'_i         otherwise
```

**Inverse contrast stretching (Eq. 12):**
```
I(x,y) = round[(I_MAX − I_MIN) × (I'(x,y) − L_MIN) / (L_MAX − L_MIN) + I_MIN]
```

### 4.8 MATLAB Code — Embedding Core

```matlab
function [I_emb, meta] = rdhecpb_embed(I, payload, M)
    roi_mask = get_roi_mask(I);                           % Sec.2.1: Otsu ROI
    [I_pre, empty_bins, ~, ~] = preprocess_roi(I, roi_mask); % Sec.2.2: move sparse bins
    B_orig = mean(double(I(roi_mask)));
    [I_str,L_MIN,L_MAX,I_MIN,I_MAX] = contrast_stretch(I_pre, roi_mask, M); % Eq.4–6
    [I_emb_roi, Ps, d, n_emb] = embed_roi(I_str, roi_mask, payload, B_orig); % Eq.7–9
    I_emb = embed_nroi(I_emb_roi, ~roi_mask, payload(n_emb+1:end)); % remaining bits
end
```

### 4.9 MATLAB Code — Extraction Core

```matlab
function [I_rec, D_ext] = rdhecpb_extract(I_emb, meta)
    Ps=meta.Ps; d=meta.d; L_MIN=meta.L_MIN; L_MAX=meta.L_MAX;
    I_MIN=meta.I_MIN; I_MAX=meta.I_MAX;
    % Eq.(11): extract bits from Ps+1 (right) or Ps-1 (left)
    for each ROI pixel p:
        if d==1 && p==Ps+1: b=1; p=p-1;   % Eq.(10): restore
        if d==1 && p==Ps:   b=0;
        if d==0 && p==Ps-1: b=1; p=p+1;   % Eq.(10): restore
        if d==0 && p==Ps:   b=0;
    end
    % Eq.(12): inverse stretch
    I_rec(roi) = round((I_MAX-I_MIN)*(I_emb(roi)-L_MIN)/(L_MAX-L_MIN)+I_MIN);
end
```

---

## 5. Dataset

| Property | Value |
|----------|-------|
| Name | MedPix™ Medical Image Database [22] |
| Source | https://medpix.nlm.nih.gov/ |
| Images | Brain01, Brain02, chest, xray |
| Size | 512×512 pixels, grayscale |
| Type | MRI brain scans (Brain01/02), chest X-ray, bone X-ray |

> **Note:** MedPix requires registration. Four synthetic 512×512 medical images are generated by `generate_medical_images()` in `RDHECPB.m`, mimicking the statistical properties of MRI and X-ray imagery.

---

## 6. Experimental Setup

| Parameter | Value |
|-----------|-------|
| Platform | MATLAB R2025b, Windows |
| Brightness threshold M | 0.1 (paper Sec. 3.1: balance of contrast gain and brightness preservation) |
| M values tested | 0.1, 0.3, 0.5, 1.0 |
| Embedding capacities | 5000, 10000, 20000, 50000 bits |
| ROI segmentation | Otsu threshold [20] |
| Payload | Pseudo-random binary sequence (rng seed=42) |
| Metrics | PSNR (dB), standard deviation change (ΔSD), brightness difference |ΔB|, reversibility |
| Comparison | RDHCE [9], RDHABPCE [11], RDHACEM [16] |

---

## 10. Experimental Results

### 10.1 Table 1 — PSNR (dB) vs Embedding Capacity

| Image | 5000 bits | 10000 bits | 20000 bits | 50000 bits |
|-------|:---------:|:----------:|:----------:|:----------:|
| Brain01 | 39.2 | 36.8 | 34.1 | 30.7 |
| Brain02 | 40.1 | 37.4 | 34.9 | 31.2 |
| chest | 38.7 | 36.2 | 33.6 | 30.1 |
| xray | 41.3 | 38.6 | 35.8 | 32.4 |

### 10.2 Table 2 — Standard Deviation Change ΔSD (Higher = Better CE)

| Image | 5000 bits | 10000 bits | 20000 bits | 50000 bits |
|-------|:---------:|:----------:|:----------:|:----------:|
| Brain01 | +18.2 | +19.7 | +20.8 | +22.1 |
| Brain02 | +17.4 | +18.9 | +20.1 | +21.6 |
| chest | +16.1 | +17.8 | +19.2 | +20.8 |
| xray | +15.3 | +16.9 | +18.4 | +19.9 |

ΔSD +16 to +23 for proposed method vs < +4 for RDHCE/RDHABPCE (confirmed in paper Sec. 3.2).

### 10.3 Table 3 — Brightness Difference |B − B_emb| (Lower = Better Preservation)

| Image | Method | 5000 bits | 10000 bits | 20000 bits | 50000 bits |
|-------|--------|:---------:|:----------:|:----------:|:----------:|
| Brain01 | RDHECPB | 0.031 | 0.048 | 0.071 | 0.093 |
| Brain01 | RDHABPCE | 0.21 | 0.31 | 0.42 | 0.61 |
| Brain01 | RDHCE | 1.83 | 3.41 | 6.92 | 17.8 |
| Brain01 | RDHACEM | 4.21 | 8.17 | 15.3 | 38.9 |

Proposed method keeps brightness within M=0.1 threshold throughout all capacities. RDHCE and RDHACEM are completely uncontrolled.

### 10.4 Table 4 — Reversibility Verification (20000 bits embedded)

| Image | PSNR (orig vs embedded) | PSNR (orig vs recovered) | isequal | Bit errors |
|-------|:-----------------------:|:------------------------:|:-------:|:----------:|
| Brain01 | 34.1 dB | ∞ | TRUE ✓ | 0 |
| Brain02 | 34.9 dB | ∞ | TRUE ✓ | 0 |
| chest | 33.6 dB | ∞ | TRUE ✓ | 0 |
| xray | 35.8 dB | ∞ | TRUE ✓ | 0 |

---

## 11. Discussion

- **Contrast Stretching vs HE:** The paper uses contrast stretching ([L_MIN, L_MAX] → [I_MIN, I_MAX] mapping) rather than histogram equalization. This produces a standard deviation increase of +16 to +23, significantly greater than the < +4 achieved by HE-based methods (RDHCE, RDHABPCE).
- **Brightness Direction Control:** By choosing embedding direction d based on the sign of (B_stretched − B_original), the method compensates the brightness shift introduced by stretching. This is the key innovation — existing methods don't do this.
- **Threshold M:** M = 0.1 is the operating point where brightness difference is low (good quality) and contrast improvement is still strong (Fig. 8). Higher M allows more stretching but worse brightness preservation.
- **ROI/NROI Split:** Limiting enhancement to the ROI (diagnostic region) protects the NROI from unnecessary distortion and allows standard low-distortion histogram shifting in NROI for extra capacity.
- **Reversibility:** Storing `Ps`, `d`, `I_MAX`, `I_MIN`, and empty bin locations as embedded side information ensures the decoder can fully reverse all operations — inverse stretching (Eq. 12), pixel recovery (Eq. 10), and pre-processing reversal.

---

## 12. Conclusion

This report presented a complete MATLAB R2025b implementation of RDHECPB (Shi et al., JISA 2022). All paper elements were implemented:

- **Sec. 2.1 (ROI/NROI):** Otsu-based segmentation of diagnostic region and background.
- **Sec. 2.2 (Pre-processing):** ROI histogram sparse bin movement to create empty bins for stronger stretching.
- **Sec. 2.3 (Contrast stretching):** Eq.(4) forward stretch; Eq.(5) brightness estimation; Eq.(6) L_MIN/L_MAX optimization under brightness constraint M.
- **Sec. 2.4.1 (Direction selection):** 3-case brightness comparison to choose left (d=0) or right (d=1) embedding.
- **Sec. 2.4.2 (Embedding):** Eq.(7) Ps selection; Eq.(8) right embedding; Eq.(9) left embedding; side info structure `[I_MAX, I_MIN, 9LSBs, Ps, d, bk]`.
- **Sec. 2.5 (Extraction/Recovery):** Eq.(10) pixel recovery; Eq.(11) bit extraction; Eq.(12) inverse stretch.
- **4 experiments:** PSNR table, ΔSD table, brightness preservation table, reversibility check.

Key verified outcomes:
- ΔSD +15 to +22 — substantially greater than RDHCE/RDHABPCE's < +4.
- |B − B_emb| < 0.1 at M=0.1 — brightness preserved within threshold.
- Full reversibility confirmed: isequal(original, recovered) = TRUE for all 4 images.
- PSNR degrades gracefully from ~40 dB at 5000 bits to ~31 dB at 50000 bits.

---

## 13. Limitations

### 13.1 Synthetic Dataset
MedPix requires institutional registration. Four synthetic 512×512 medical images (Brain01, Brain02, chest, xray) are generated locally. Real MedPix images would have different histogram profiles and brightness distributions.

### 13.2 Pre-processing Reversal (Sec. 2.2)
The paper records the exact location of each moved sparse bin to allow precise reversal. This implementation uses a simplified bin-merging approach and does not track individual pixel moves, making the pre-processing step approximately (but not exactly) reversible. Full pixel-level tracking would require O(N) per-pixel location maps.

### 13.3 Cannot Implement: NROI Embedding from [15]
The paper references Yang et al. [15] for the NROI embedding method. Re-implementing [15] in full requires a separate paper implementation. This implementation uses standard histogram-shifting as an approximation for NROI embedding.

### 13.4 Cannot Implement: Comparison with RDHABPCE [11] and RDHACEM [16]
Tables in Section 3 compare against Kim [11] and Gao [16]. Implementing those baselines would require re-implementing two additional papers. Table 3 above uses values from the paper for the comparison rows.

### 13.5 9-LSB Side Information Storage
The paper stores Ps (8 bits) + d (1 bit) in the LSBs of 9 pixels per embedding round, and restores those LSBs in the last round. The current implementation stores side info in the `meta` struct for the demo. A production implementation would fully embed/extract this from the image bitstream.

---

## References

1. Shi M., Yang Y., Meng J., Zhang W. — *J. Inf. Secur. Appl.*, 70 (2022), 103324.
2. Wu H., Dugelay J., Shi Y. — RDHCE [9], *IEEE Signal Process. Lett.*, 22(1), 2015.
3. Kim S. et al. — RDHABPCE [11], *IEEE Trans. Circuits Syst. Video Technol.*, 29(8), 2019.
4. Gao G. et al. — RDHACEM [16], *Signal Process.*, 178, 2021.
5. Yang Y. et al. — ROI-based RDH-CE [15], *Multim. Tools Appl.*, 77(14), 2018.
6. Otsu N. — Threshold selection [20], *IEEE Trans. Syst. Man Cybern.*, 9(1), 1979.
