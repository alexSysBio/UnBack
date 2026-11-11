# UnBack

This repository includes methods for subtracting the background in fluorescence microscopy images of cells
(e.g. bacteria on agarose pads). The background is estimated from the cell-free pixels only:
- The image is divided into square tiles (sectors).
- In each tile, the median of the pixels away from any cell is taken.
- These medians are smoothed into a background surface, which is subtracted from the image.

This tile implementation is fast. A rolling window with a 1-pixel step would be an alternative, but a slower one.

Two versions are included:

| file | version | |
|---|---|---|
| `unback2.py` | **2 (recommended)** | same method and interface, with the fixes listed in [What changed in version 2](#what-changed-in-version-2) |
| `uneven_background_correction.py` | 1 | the original implementation, kept unchanged so earlier analyses can be reproduced |

Other files:

| file | |
|---|---|
| `nd2_to_array.py` | reads Nikon `.nd2` files into numpy arrays (single images, time-lapses, multiple positions and channels) |
| `apply_background_correction_nisjobs.py` | batch example: background of every frame of every position of a NIS-Elements JOBS acquisition |
| `test_uneven_background_correction.ipynb` | step-by-step example of version 1 |

---

## Installation

Clone or download this repository and add its folder to the Python path.

Requirements:
- Python 3, `numpy`, `scipy`, `scikit-image`;
- `matplotlib` for `show=True`;
- `pims` and `pims_nd2` for `nd2_to_array.py`.

Tested with numpy 2.1, scipy 1.15 and scikit-image 0.25.

```python
import sys
sys.path.append(r'path/to/UnBack')
import unback2 as bkg
```

---

## How the background is estimated

1. **Cells.** The cells are found in the phase-contrast image: `get_inverted_mask()` inverts it (cells are dark
   in phase) and applies an Otsu threshold. In version 2, segmentation masks (e.g. Omnipose or Cellpose
   labels) can be added with `cell_masks`.
2. **Margin.** Pixels within `dilation` px of a cell are left out. Bright cells are surrounded by a glow of
   out-of-focus and scattered light, which must not be counted as background. The pixels that remain are
   the cell-free pixels.
3. **Tiles.** The image is divided into square tiles of `estimation_step` px, and each tile gets the median of
   its cell-free pixels.
4. **Surface.** The tile medians are smoothed with a Gaussian (`smoothing_sigma` px) into a background surface.
5. **Subtraction.** `subtracted = image − surface`.
6. **Flat field.** The surface also follows the excitation light, because the background fluorescence
   (medium, agarose, glass) is excited by the same light as the cells. Dividing by the normalised surface
   corrects uneven excitation: `corrected = subtracted / profile`.

---

## How to use

### 1. Read the images

```python
import numpy as np
import tifffile
import nd2_to_array as ndtwo
import unback2 as bkg

# TIF files
phase_image = tifffile.imread('Point0000_ChannelTrans.tif').astype(float)
fluor_image = tifffile.imread('Point0000_ChannelmCherry.tif').astype(float)

# or ND2 files (single image per file: element [2] is the image)
phase_image = ndtwo.nd2_to_array('Point0000_ChannelTrans_Seq0000.nd2')[2].astype(float)
fluor_image = ndtwo.nd2_to_array('Point0000_ChannelmCherry_Seq0001.nd2')[2].astype(float)
```

What element `[2]` of `nd2_to_array` holds depends on the acquisition:

| acquisition | `[2]` |
|---|---|
| single image | the image |
| stream / time-lapse (`'t'`) | `{frame: image}` |
| several positions and channels (`'mc'`) | `{position: {channel: image}}` |

### 2. Make the cell mask

```python
phase_mask = bkg.get_inverted_mask(phase_image)      # True = cell (dark core in phase contrast)
```

If you have segmentation masks, pass them too. The Otsu mask finds only the dark core of each cell, about
65 % of its pixels; the segmentation masks cover whole cells. The Otsu mask in turn also finds cells the
segmentation missed.

```python
labels = tifffile.imread('Point0000_ChannelTrans_cp_masks.tif')   # e.g. Omnipose / Cellpose labels
```

Check the mask by eye:

```python
import matplotlib.pyplot as plt
plt.imshow(phase_image, cmap='gray'); plt.contour(phase_mask, [0.5], colors='r', linewidths=0.5); plt.show()
```

### 3. Subtract the background

```python
corrected, subtracted, surface, excluded = bkg.back_sub(
    fluor_image, phase_mask,
    dilation=30,              # px around every cell left out of the background
    estimation_step=128,      # px, side of the tiles
    smoothing_sigma=100,      # px, smoothing of the tile medians
    cell_masks=labels)        # optional
```

| output | |
|---|---|
| `corrected` `[0]` | background subtracted and divided by the flat-field profile (uneven excitation corrected) |
| `subtracted` `[1]` | background subtracted |
| `surface` `[2]` | the background surface (counts) |
| `excluded` `[3]` | pixels left out of the background (True = cell, or within `dilation` px of one) |

Which output to use:
- `subtracted` when only the background has to go: segmentation, spot detection, intensities within a cell,
  ratios between compartments of the same cell.
- `corrected` when absolute intensities are compared between cells in different parts of the field (e.g.
  concentrations). This needs the right `camera_offset`; see [step 6](#6-flat-field-and-camera-offset).

### 4. Check the result

```python
out = bkg.back_sub(fluor_image, phase_mask, cell_masks=labels, return_info=True, show=True)
corrected, subtracted, surface, excluded, info = out

residual = bkg.background_pixels(subtracted, excluded)
print(np.median(residual))        # ~0 if the surface is right
print(info['noise_sd'])           # pixel noise SD (e.g. for signal-to-noise ratios)
print(info['n_empty_tiles'])      # tiles without enough cell-free pixels, filled in from their neighbours
print(info['free_fraction'])      # fraction of the image used as background
```

`show=True` plots the cell-free pixels, the surface and the corrected image.

### 5. Choose the settings

**`dilation`: cover the glow of the cells.** Measure how far the glow reaches:

```python
from scipy import ndimage
cells = phase_mask | (labels > 0)
dist = ndimage.distance_transform_edt(~cells)              # px to the nearest cell
far = np.median(fluor_image[dist > 60])
for a, b in [(1, 3), (3, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 40), (40, 60)]:
    ring = (dist > a) & (dist <= b)
    print(f'{a:>2}-{b:<2} px: {np.median(fluor_image[ring]) - far:7.1f} counts above the far background')
```

- The glow fades slowly, so `dilation` is a compromise: a wider margin removes more glow but leaves fewer
  cell-free pixels.
- Example: bright HupA-mCherry cells at 66 nm/px, with 'far' meaning pixels more than 60 px from any cell.
  The excess was 290 counts at 5–10 px, 190 at 10–15 px, 80 at 20–30 px, 44 at 30–40 px and 18 at 40–60 px.
  - Widening the margin from 30 to 60 px lowered the background by 8–15 counts and changed cell intensities
    by only +0.3 to +0.6 %.
  - But it left only 25–46 % of the field as background, and 3–27 of the 256 tiles were empty.
  - 30 px (63–77 % of the field cell-free, no empty tiles) was the better compromise.
- The glow scales with the brightness of the cells, so the relative error is similar for dim cells, and dim
  signals need a smaller margin.
- A margin that is too large in a crowded field leaves tiles without cell-free pixels. These are filled in
  from their neighbours (`info['n_empty_tiles']`). If no tile has `min_pixels` left, `back_sub` raises an
  error: use a smaller `dilation`.

**`estimation_step` and `smoothing_sigma`: how fast the background may change across the field.**
- Large values follow only slow changes, such as uneven illumination, and are robust against local
  contamination.
- Small values follow faster changes, but also the glow of dense groups of cells.
- The defaults (128 px tiles, σ = 100 px) suit a 2048 × 2048 field.
- The tile size doesn't have to divide the image size.

**`min_pixels`** (default 50): cell-free pixels a tile needs for its median to count.

### 6. Flat field and camera offset

Dividing by the surface corrects uneven excitation only if two things hold:
- The background fluorescence scales with the excitation light.
- The camera offset is removed first. Most scientific cameras add a fixed offset, often 100 counts, to every
  pixel. It doesn't scale with the light, so it must not be part of the excitation profile.

Measure the offset as the median of an image taken without light (shutter closed, same exposure and binning),
or take it from the camera specification:

```python
camera_offset = np.median(dark_image)
corrected = bkg.back_sub(fluor_image, phase_mask, cell_masks=labels, camera_offset=camera_offset)[0]
```

- The profile is `(surface − camera_offset)` divided by its median. Corrected intensities therefore stay on
  the scale of the raw image and can be compared between images.
- With `camera_offset=0`, the profile is flatter than the real excitation and the correction is too weak.
- In a synthetic test (excitation 30 % lower in the corners, offset 100, identical cells), the spread of the
  corrected cell intensities (IQR / median) was 0.048 with `camera_offset=0` and 0.012 with
  `camera_offset=100`.

To check whether the flat field helps in your images: cells sitting on a higher background should be brighter,
and `corrected` should make the cells more alike than `subtracted` does. This is a weak test when cells differ a
lot among themselves.

```python
from scipy.stats import spearmanr
ids = np.unique(labels)[1:]
m_sub = ndimage.mean(subtracted, labels, ids)
m_cor = ndimage.mean(corrected, labels, ids)
s = ndimage.mean(surface, labels, ids)
spread = lambda v: np.subtract(*np.percentile(v, [75, 25])) / np.median(v)
print(spearmanr(m_sub, s)[0], spread(m_sub), spread(m_cor))   # rho > 0 and spread(m_cor) < spread(m_sub)
```

### 7. Time-lapses and streams

The cells barely move during a stream, so the mask is made once and every frame is corrected:

```python
frames = ndtwo.nd2_to_array('Point0000_ChannelFITC_Seq0002.nd2')[2]       # {frame: image}
phase_mask = bkg.get_inverted_mask(phase_image)

subtracted_frames, noise = {}, {}
for t, frame in frames.items():
    out = bkg.back_sub(frame.astype(float), phase_mask, cell_masks=labels, return_info=True)
    subtracted_frames[t] = out[1]
    noise[t] = out[4]['noise_sd']
```

This takes about 2.7 s per 2048 × 2048 frame. For long streams where the background changes only in level
(e.g. bleaching), not in shape, the surface can be estimated once on the mean image and shifted per frame:

```python
mean_image = np.mean([f.astype(float) for f in frames.values()], axis=0)
_, _, surface, excluded = bkg.back_sub(mean_image, phase_mask, cell_masks=labels)
free = ~excluded
ref = np.median(mean_image[free])
subtracted_frames = {t: f - (surface + np.median(f[free]) - ref) for t, f in frames.items()}
```

### 8. Several positions and channels in one ND2

```python
axis, _, images, channels, *_ = ndtwo.nd2_to_array('snapshots.nd2')       # axis 'mc'
for pos in images:                                                       # {position: {channel: image}}
    phase_mask = bkg.get_inverted_mask(images[pos]['Trans'])
    out = bkg.back_sub(images[pos]['mCherry'].astype(float), phase_mask)
```

Channel names (`'Trans'`, `'mCherry'`, …) are those of the acquisition; see `channels`.

### 9. Switching from version 1

```python
import unback2 as bkg          # instead of: import uneven_background_correction as bkg
```

- **Same calls.** `back_sub(signal_image, phase_mask, dilation, estimation_step, smoothing_sigma, show)` keeps
  its arguments and their order, so positional calls work unchanged. This includes
  `apply_background_correction_nisjobs.py`, which passes `dilation=15` explicitly. Consider 30 for bright
  signals.
- **The `[0]` outputs change scale.** Version 2 normalises the flat field by the median instead of the
  maximum. Values from versions 1 and 2 shouldn't be mixed: re-run all images when switching.
- **The `[1]` outputs barely change.** Subtraction alone changes by less than 1 % for bright cells.
- **Reproducing version 1.** `flat_norm='max', camera_offset=0` reproduces version 1, apart from the shape of
  the margin. With the same excluded pixels, the result is identical.
- **Noise.** For the background noise, use `info['noise_sd']` (`return_info=True`) instead of the SD of the
  background pixels.

---

## Reference

### `unback2.py`

`back_sub(signal_image, phase_mask=None, dilation=30, estimation_step=128, smoothing_sigma=100, show=False,
cell_masks=None, camera_offset=0.0, flat_norm='median', min_pixels=50, return_info=False)`

| parameter | default | meaning |
|---|---|---|
| `signal_image` | – | the fluorescence image (2-D) |
| `phase_mask` | None | boolean cell mask, e.g. `get_inverted_mask(phase_image)` |
| `cell_masks` | None | more cells to exclude: label image or boolean mask; at least one of `phase_mask` / `cell_masks` is needed |
| `dilation` | 30 | px (Euclidean distance) around the cells left out of the background |
| `estimation_step` | 128 | px, side of the square tiles |
| `smoothing_sigma` | 100 | px, Gaussian smoothing of the tile medians (0: none, piecewise constant) |
| `min_pixels` | 50 | cell-free pixels a tile needs to count |
| `camera_offset` | 0 | counts the camera adds without light; removed from the flat-field profile only |
| `flat_norm` | `'median'` | normalisation of the flat-field profile: `'median'`, or `'max'` (version 1) |
| `show` | False | plot the cell-free pixels, the surface and the corrected image |
| `return_info` | False | also return a dict of diagnostics |

Returns `(corrected, subtracted, surface, excluded)`, plus `info` with `return_info=True`. The keys of `info`:
- `noise_sd`: pixel noise SD, from differences of neighbouring cell-free pixels.
- `background_median`: median of the surface.
- `surface_min`, `surface_max`.
- `flat_min`, `flat_max`: range of the flat-field profile.
- `free_fraction`: fraction of the image that is cell-free.
- `n_tiles`, `n_empty_tiles`.
- `dilation`.

| function | |
|---|---|
| `get_inverted_mask(phase_image)` | Otsu mask of the inverted phase image (True = cell) |
| `cell_free_pixels(cells, dilation)` | True where a pixel is more than `dilation` px from every cell |
| `tile_medians(image, free, step, min_pixels)` | median and number of cell-free pixels of every tile |
| `cell_free_bkg_estimation(image, free, step, sigma, min_pixels)` | the background surface, plus the tile medians and counts |
| `noise_sd(image, free)` | pixel noise SD: 1.4826 × MAD of neighbouring-pixel differences / √2 |
| `background_pixels(corrected, excluded)` | the cell-free pixels of a corrected image |

### `uneven_background_correction.py` (version 1)

`back_sub(signal_image, phase_mask, dilation=15, estimation_step=128, smoothing_sigma=60, show=False)` returns
`(corrected, subtracted, surface, excluded)`:
- the cell mask is grown with `binary_dilation(iterations=dilation)`;
- the flat-field profile is the surface divided by its maximum.

`get_inverted_mask(phase_image)` and `cell_free_bkg_estimation(masked_signal_image, step)` work as described in
the code.

---

## What changed in version 2

Version 2 (2026) keeps the method and the interface of version 1 (2024).

**Circular margin around cells.** Cells are excluded by Euclidean distance.
- `binary_dilation(iterations=n)` gave a diamond that reaches only n/√2 px diagonally (21 px for 30), so the
  glow of bright cells leaked into the background where cells are dense.
- On 6 fields of HupA-mCherry cells, the glow left in the surface fell from 6–12 to 3–6 counts (median) and
  from 16–81 to 9–29 counts (maximum).

**Segmentation masks.** `cell_masks` (e.g. Omnipose labels) excludes whole cells together with the Otsu mask of
the phase image, which finds only the dark core of each cell (about 65 % of its pixels). Either mask can be
used alone.

**No more NaN images.** A tile without enough cell-free pixels is left out and filled in by normalised
Gaussian smoothing.
- Before, its NaN median spread over the image: one empty 128 px tile made 32 % of a 512 × 512 image NaN.
- When every tile has pixels, the surface is identical to version 1.

**Flat field normalised by the median.** The flat-field profile is normalised by its median instead of its
maximum.
- With the maximum, each image was scaled to its brightest background spot, so corrected intensities
  differed between images of the same sample.
- Between the 6 positions of one experiment, the spread was 20 %; it is now 1 %.

**Camera offset.** `camera_offset` removes the camera's dark offset from the flat-field profile, because the
offset doesn't scale with the excitation light. In a synthetic test this made the corrected cell intensities
four times more uniform (spread 0.048 → 0.012).

**Noise estimate.** `return_info=True` returns the noise SD, from differences of neighbouring cell-free pixels.
- This is blind to glow and uneven background: 17–18 counts, against 29–38 for the SD of the background
  pixels.
- It also returns diagnostics: empty tiles, cell-free fraction, surface and flat-field range.

**Boolean masks throughout.** Cell-free pixels with value 0 are no longer dropped from the background.

**Tile size.** The tiles no longer have to divide the image size.

**Backward compatibility.** `flat_norm='max'` and `camera_offset=0` reproduce version 1, apart from the shape of
the margin; with the same excluded pixels the result is identical. Speed is about the same (2.7 s against
2.5 s per 2048 × 2048 image).

---

## Limitations

- Glow from a neighbouring cell that lands *on* a cell isn't background and isn't removed. Intensities of cells
  in tight groups include some of their neighbours' light.
- The surface follows only changes slower than about `smoothing_sigma`.
- `get_inverted_mask` assumes cells are darker than the background in phase contrast, and that the phase image
  is even enough for a single Otsu threshold. Otherwise, pass segmentation masks as `cell_masks`.
- The flat-field correction assumes the background fluorescence scales with the excitation light (see
  [step 6](#6-flat-field-and-camera-offset)).

---

## Cite

Nonequilibrium polysome dynamics promote chromosome segregation and its coupling to cell growth in
*Escherichia coli*.

Alexandros Papagiannakis, Qiwei Yu, Sander K. Govers, Wei-Hsiang Lin, Ned S. Wingreen, Christine Jacobs-Wagner.

eLife 2025;14:RP104276. DOI: https://doi.org/10.7554/eLife.104276.3

## License

MIT (see `LICENSE`).
