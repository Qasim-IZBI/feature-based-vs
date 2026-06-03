# Research Discussion: Registration-Informed Conditioning for Virtual Staining Diffusion Models

## Context

### MIU-Diff (Score-based Diffusion Model for Unpaired Virtual Histology Staining)

**Task:** H&E → IHC virtual staining (unpaired)

**Architecture:** Two-stage score-based diffusion model
- **Stage 1:** Unconditional diffusion model pretrained on IHC images alone — learns marginal distribution `q(y₀)`
- **Stage 2:** MI-guided reverse diffusion conditioned on the H&E input — learns conditional distribution `q(y₀|x₀)`

**Current conditioning signal:**
- Color is stripped from H&E: `x₀ → x'₀` (grayscale)
- A **gradient map** `g_{x,t}` (Sobel edges) of the grayscale H&E serves as the structural proxy
- The MI estimator `G_θ` is trained entirely within the IHC domain: it measures shared information between IHC image `y` and its own gradient map `g_y`
- During inference, the H&E gradient map substitutes for the IHC gradient map

**Three key components:**
1. Global MI-guided energy function `M(y, x, t)` — disentangles staining style from tissue structure
2. Timestep-adaptive reverse-time SDE — controllable staining intensity and structural reconstruction
3. Local MI-driven contrastive learning `ℓ_PCL` — patch-level structural consistency

**Limitation of current conditioning:**
The gradient map captures *all* edges in H&E — including staining-specific color transitions, blood artifacts, and noise — not just structures genuinely shared with IHC. It is a proxy for shared content, not a direct measurement of it.

---

## Registration Papers Overview

### 1. VALIS — Virtual Alignment of pathoLogy Image Series (Gatenbee et al., Nature Comms 2023)
- Registers multiplex WSI (H&E, IHC, IF) — any number of slides
- Pipeline: Convert → Normalize → Feature detection → Match → Sort → Serial rigid → Non-rigid → Micro-registration
- Uses both deep learning features (VGG) and handcrafted (BRISK) descriptors
- Maximizes mutual information between matched features at convergence
- Can register different modalities including H&E and IHC

### 2. DeeperHistReg (Wodzinski et al., 2024)
- Three modules: preprocessing, initial alignment (affine), deformable registration
- Preprocessing: grayscale conversion + CLAHE normalization
- Initial alignment: SuperPoint + SuperGlue (self-supervised, no stain-specific training)
- Deformable: B-Splines or dense displacement fields

### 3. RegWSI — Winner of ACROBAT 2023 (Wodzinski et al., Computer Methods 2024)
- Two-step hybrid: (i) SuperPoint/SuperGlue initial alignment + (ii) intensity-based nonrigid registration
- Similarity measure: Local Normalized Cross-Correlation (NCC)
- Deformation model: B-Splines
- No fine-tuning required — generalizes across tissue types and stains
- Incorporated into DeeperHistReg framework

### 4. STAR — Serial Tissue Alignment for Rigid Registration (Liu & Ding, 2025)
- Fast, lightweight rigid registration framework
- **Stain-conditioned preprocessing** — key insight: H&E gets histogram equalization + intensity inversion; IHC gets background thresholding + Gaussian blur + intensity inversion → both modalities mapped to similar intensity space
- Multi-stage coarse-to-fine template matching via 2D cross-correlation
- Works across H&E, IHC, PAS, PASM, special stains

### 5. CORE — Cell-Level Coarse-to-Fine Registration (Nasir et al., 2025)
- Two-stage framework:
  - **Coarse:** prompt-based tissue mask → dense feature matching → rigid + non-rigid alignment
  - **Fine:** nuclei detection on both modalities → shape-aware point-set registration → CPD (Coherent Point Drift) for non-rigid displacement field
- Key insight: **nuclei centroids** are the most intrinsic, stain-agnostic features across H&E, IHC, PAS, and mIF
- Achieves nuclei-level correspondence across modalities

### 6. UWarp (Schieb et al., 2025)
- Designed for same-stain, different-scanner registration
- Hierarchical: global affine + fine local non-rigid corrections
- Automatic landmark generation and matching with quality scoring
- Patch-level alignment accuracy — useful for patch-level domain shift analysis

---

## The Core Idea

### Problem Statement

Replace or augment MIU-Diff's gradient map condition with a **registration-informed shared-structure representation** computed from H&E alone — one that encodes what registration methods show is genuinely cross-modal invariant, rather than all edges indiscriminately.

### Key Insight from Registration Literature

Across all registration papers, the most reliable cross-modal anchors are anatomical structures that are intrinsically stain-agnostic:

| Registration Method | Shared Feature Exploited |
|---|---|
| VALIS | VGG + BRISK keypoints on normalized grayscale |
| DeeperHistReg / RegWSI | SuperPoint keypoints (self-supervised) + NCC |
| STAR | Stain-conditioned preprocessing mapping both modalities to similar intensity space |
| UWarp | Automatic landmark generation at tissue boundaries |
| CORE | **Nuclei centroids** — detected independently in each modality |

---

## Proposed Conditioning Signals

### Stain-Dependent Strategy

The right shared-structure conditioning depends on the target stain:

| Target Stain | What is shared with H&E | Appropriate conditioning |
|---|---|---|
| IHC | Nuclei (hematoxylin counterstain present) | Nuclei feature map |
| PSR (Picrosirius Red) | Collagen/stromal architecture — **NOT nuclei** | Tissue compartment map, stromal regions |
| PAS | Gland/mucosa boundaries | Tissue architecture map |
| mIF | Channel-specific cell populations | Cell boundary / compartment map |

### Important Stain-Specific Note — PSR

**Picrosirius Red (PSR) does not highlight nuclei.** PSR targets collagen fibers — nuclei are essentially invisible. Therefore nuclei-based conditioning breaks down for H&E → PSR translation. The shared content between H&E and PSR is:
- Collagen-dense stromal regions
- Tissue compartment boundaries (epithelial vs. stromal interfaces)
- Vessel walls and gland boundaries

---

## Conditioning Signal Options (No Labels, No AI Model Required)

### Option A — Classical Thresholding / Algebraic Methods

**Color Deconvolution (Ruifrok & Johnston, 2001)**
- Mathematically separates H&E into hematoxylin and eosin channels using known optical density matrix
- Hematoxylin channel → nuclei / nuclear-dense regions (epithelium)
- Eosin channel → cytoplasm, ECM, stroma
- Otsu threshold on each channel gives rough compartment masks
- Available in scikit-image, histomicstk

**K-means Clustering in LAB/HSV Color Space**
- H&E pixels cluster naturally into biologically meaningful groups
- K=3: background/lumen, nuclear/epithelial, stroma/cytoplasm
- K=4-5: further separates lumen from background, stroma from cytoplasm
- Completely unsupervised, no labels, no model

**HSV Direct Thresholding**

| Compartment | Hue Range (approx.) |
|---|---|
| Nuclei | 200–280° (blue-purple) |
| Eosin / stroma | 300–360° (pink-red) |
| Background / lumen | Any, high Value (>0.9) |

Stable across slides after basic stain normalization (Macenko/Vahadane — both classical SVD-based, no AI).

**For PSR specifically:** Simple threshold on red channel saturation in HSV reliably separates collagen-dense stroma (high red saturation) from cellular regions and background.

### Option B — Registration-Derived Classical Representations

These are more powerful than thresholding because they are specifically designed to be **cross-modal invariant** — they suppress staining-specific intensity information and retain only structural information.

**Normalized Gradient Field (NGF)**

Used in HistokatFusion (ANHIR winner). Formula:
```
NGF = ∇I / (|∇I| + ε)
```
- Divides gradient by its magnitude → removes intensity information, keeps only **orientation**
- Robust to intensity differences between H&E and IHC/PSR
- Purely algebraic, no model needed
- Direct improvement over the raw Sobel gradient currently used in MIU-Diff

**Difference of Gaussians (DoG) — multi-scale**
- Underlying detection mechanism of SIFT
- Highlights structurally stable locations across scale: gland boundaries, vessel walls, tissue interfaces
- Intensity-independent by construction
- Classical computation, no AI

**Phase Congruency (Kovesi, 1999)**
- Detects features where Fourier components are in phase
- Entirely **intensity-invariant** — responds to structural transitions regardless of staining intensity
- Used in cross-modal medical image registration for this exact property
- More computationally intensive but classical (FFT-based)

**Histogram of Oriented Gradients (HOG)**
- Dense descriptor from local gradient orientations
- Captures tissue texture and boundary structure
- Relatively insensitive to absolute intensity values

### Comparison of Conditioning Signals

| Signal | Intensity-Invariant | Cross-Modal Stable | Complexity |
|---|---|---|---|
| Sobel gradient (current MIU-Diff) | No | No | Trivial |
| Color deconvolution | N/A | Domain-specific | Algebraic |
| K-means clusters | Partial | Partial | Algebraic |
| NGF | Yes | Yes | Algebraic |
| DoG (multi-scale) | Partial | Yes | Classical |
| Phase congruency | Yes | Yes | Classical (FFT) |
| HOG | Partial | Yes | Classical |

---

## Proposed Integration into MIU-Diff

### Minimal Change (most testable)

Replace the plain Sobel gradient map with **NGF**:
```
condition: g_{x,t} (Sobel)  →  NGF_{x,t}
```
The MI estimator `G_θ` is retrained using `NGF(y)` instead of `g_y` on the IHC side. The unpaired training framework is fully preserved. This directly tests whether a cross-modal-invariant gradient representation is a better condition than a raw gradient map.

### Extended Conditioning

Stack multiple channels as the condition:
```
condition = [NGF(x'₀), DoG(x'₀), tissue_compartment_map(x'₀)]
```
Where `tissue_compartment_map` is computed via color deconvolution + thresholding.

### Registration-Supervised Training Signal

During training, run a registration method (STAR for speed, CORE for nuclei-level accuracy) on H&E-IHC pairs to find matched landmark positions. Use these to:
1. Construct a **landmark confidence map** — spatial heatmap over H&E indicating where structural correspondence to IHC is strongest
2. Re-weight the contrastive loss `ℓ_PCL` — upweight patches centered on high-confidence landmarks

### For IHC Target — Nuclei-Based Conditioning

1. Detect nuclei in H&E using classical watershed on the hematoxylin channel (after color deconvolution)
2. Generate nuclei centroid map + local morphological features (area, density, eccentricity)
3. Use as conditioning channel — directly equivalent to what CORE uses as fine-registration landmarks
4. MI estimator trained on IHC side using nuclei detected from IHC hematoxylin counterstain

---

## Open Questions / Next Steps

- Which combination of classical conditioning signals gives the best MI estimate with IHC?
- Does NGF alone outperform Sobel gradient in the ablation, or is the multi-channel stack needed?
- For PSR: does a stromal probability map (from color deconvolution) provide useful conditioning?
- How to handle stains where neither nuclei nor clear compartment boundaries exist?
- Can landmark confidence maps from registration be practically computed at training scale for WSI patches?
