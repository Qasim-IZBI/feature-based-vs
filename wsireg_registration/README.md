# H&E ↔ PSR WSI Registration — Mouse Liver

Registers H&E and PSR (Picrosirius Red) stained whole slide images from the
same mouse liver sample. Sections may be 2–3 cross-sections apart.

## Pipeline

```
rigid → affine → non-linear (B-spline)
```

| Stage | What it corrects |
|---|---|
| Rigid | Slide placement — rotation + translation |
| Affine | Cutting angle — shear + global scale |
| Non-linear | Local tissue deformation from section gap |

**Similarity metric:** Mutual Information (multimodal-safe across H&E ↔ PSR)  
**Fixed / reference:** PSR  
**Moving:** H&E warped into PSR space

## Setup

### 1. Install the OpenSlide C library

```bash
# macOS
brew install openslide

# Ubuntu / Debian
sudo apt install libopenslide-dev
```

### 2. Install Python dependencies

**Option A — pip (Intel Mac / Linux)**
```bash
pip install -r requirements.txt
```

**Option B — conda (recommended on Apple Silicon)**
```bash
# itk-elastix has no ARM wheel on PyPI; conda-forge has one
conda install -c conda-forge itk-elastix
pip install wsireg openslide-python tifffile matplotlib Pillow
```

## Usage

```bash
# Basic registration
python register.py --he path/to/liver_HE.svs --psr path/to/liver_PSR.svs

# With QC overlay PNG
python register.py --he liver_HE.svs --psr liver_PSR.svs --qc

# Custom output directory and pixel size (20× scanner = 0.46 µm/px)
python register.py --he liver_HE.svs --psr liver_PSR.svs \
    --out ./results --pixel-size 0.46 --qc
```

Output files in `./registration_output/` (or `--out`):
- `*HE_registered_to_PSR*.ome.tiff` — registered H&E as a pyramidal OME-TIFF
- `qc_overlay.png` — four-panel QC: PSR | H&E original | H&E registered | checkerboard

The OME-TIFF can be opened alongside the original PSR in QuPath, OMERO, or Fiji/BigDataViewer.

## Key parameters to tune

### Pixel size (`--pixel-size`)
Check your scanner's metadata — wrong pixel size will cause Elastix to use the
wrong physical grid spacing for the non-linear step.

| Scanner / objective | Typical µm/px |
|---|---|
| Leica Aperio 40× | 0.23 |
| Leica Aperio 20× | 0.46 |
| Hamamatsu NDP 40× | 0.23 |
| Zeiss Axio Scan 20× | 0.44 |

### Preprocessing channels (`HE_PREPRO` / `PSR_PREPRO` in `register.py`)

Both stains use the **green channel** by default. This works well because:
- H&E: nuclei and cytoplasm absorb green → dark tissue on bright background
- PSR: red collagen absorbs green → dark collagen; hematoxylin counterstain shows nuclei

If registration quality is poor (e.g. very faint counterstain in PSR), try:
```python
# In register.py — switch to grayscale average across all channels
HE_PREPRO  = {"as_uint8": True, "ch_indices": [0, 1, 2]}
PSR_PREPRO = {"as_uint8": True, "ch_indices": [0, 1, 2]}
```

### Non-linear grid spacing
If local deformation from 3+ section gaps is too large for the default B-spline
grid, add a custom Elastix parameter map:

```python
from wsireg.parameter_maps.reg_params import get_elastix_parameter_map

nl_params = get_elastix_parameter_map("nl")
nl_params["FinalGridSpacingInPhysicalUnits"] = ["50.0"]  # µm; default ~100
reg.add_reg_path("HE", "PSR", reg_params=["rigid", "affine", nl_params])
```

## If registration quality is insufficient

For sections more than 3 apart with severe warping, consider
**DeepHistReg** (DL-based, handles larger deformations):
```bash
pip install deephistreg
```
It can be used as a drop-in replacement for the non-linear stage.
