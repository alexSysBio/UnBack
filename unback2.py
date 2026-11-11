# -*- coding: utf-8 -*-
"""
UnBack2: cell-free background estimation and subtraction for fluorescence images.

Version 2 of uneven_background_correction.py (UnBack). The method is the same -- the median of the cell-free
pixels in square tiles, smoothed with a Gaussian into a background surface that is subtracted -- and back_sub()
takes the same arguments and returns the same four images, so existing code can switch with

    import unback2 as bkg          # instead of: import uneven_background_correction as bkg

Changes from version 1:

1. Circular margin. Pixels closer than `dilation` px (Euclidean distance) to a cell are excluded. Version 1
   grew the cell mask with binary_dilation(iterations=dilation), which gives a diamond whose diagonal reach
   is only dilation / sqrt(2) (21 px for 30), while the glow of bright cells (out-of-focus and scattered
   light) still adds ~60 counts at 20-30 px (HupA-mCherry, stationary phase).
2. Cell masks. `cell_masks` (e.g. the Omnipose labels) are excluded together with the Otsu mask of the
   inverted phase image (`phase_mask`); either can also be given alone.
3. Empty tiles. Tiles with fewer than `min_pixels` cell-free pixels are left out, and the surface is computed
   by normalised Gaussian smoothing (smoothed tile values / smoothed tile weights), so the other tiles fill
   them in. In version 1 such a tile gives a NaN median, and the smoothing spreads the NaN over a large part
   of the image (one empty 128 px tile in a 512 px image: 32 % of the image NaN). When every tile has enough
   pixels the surface is the same as in version 1.
4. Flat field. The flat-field corrected image [0] = (signal - surface) / profile, with the profile
   (surface - camera_offset) normalised by its MEDIAN. Version 1 normalised by the maximum, so every image was
   scaled to its single brightest background spot and corrected intensities differed between images of the
   same sample (8-29 % between the positions of one experiment). The camera offset (counts the camera adds
   without light; 0 by default) does not scale with the excitation light and is removed from the profile
   first. flat_norm='max' with camera_offset=0 reproduces version 1.
5. Noise. back_sub(..., return_info=True) also returns the noise SD of the background, from differences of
   neighbouring cell-free pixels (1.4826 x MAD / sqrt 2), which is blind to the glow of the cells and to the
   uneven background (the SD of the background pixels is not: 29-38 counts against a pixel noise of 18-19).
6. Boolean masks throughout: cell-free pixels whose value is 0 are no longer dropped.

Based on uneven_background_correction.py (Alexandros Papagiannakis, HHMI at Stanford University, 2024).
"""

import numpy as np
from scipy import ndimage
from skimage.filters import threshold_otsu

DILATION = 30            # px (Euclidean) around every cell excluded from the background
ESTIMATION_STEP = 128    # px, side of the square tiles
SMOOTHING_SIGMA = 100    # px, Gaussian smoothing of the tile medians
MIN_PIXELS = 50          # cell-free pixels a tile needs to count


# ---------------------------------------------------------------- cell masks
def get_inverted_mask(phase_image):
    """Cells of a phase-contrast image: Otsu threshold of the inverted image (cells are dark).
    -> boolean 2D array, True = cell"""
    phase = np.asarray(phase_image, float)
    positive = phase[phase > 0]
    floor = positive.min() if positive.size else 1.0
    inverted = 1.0 / np.where(phase > 0, phase, floor)
    return inverted > threshold_otsu(inverted.ravel())


def cell_free_pixels(cells, dilation=DILATION):
    """True where a pixel is more than `dilation` px (Euclidean distance) from every cell pixel."""
    cells = np.asarray(cells, bool)
    if not cells.any():
        return np.ones(cells.shape, bool)
    if dilation <= 0:
        return ~cells
    return ndimage.distance_transform_edt(~cells) > dilation


# ---------------------------------------------------------------- background surface
def tile_medians(signal_image, free, step=ESTIMATION_STEP, min_pixels=MIN_PIXELS):
    """Median of the cell-free pixels of every step x step tile (NaN if fewer than min_pixels).
    -> (medians, pixel counts), arrays of tiles"""
    H, W = signal_image.shape
    ny, nx = -(-H // step), -(-W // step)
    med = np.full((ny, nx), np.nan)
    count = np.zeros((ny, nx), int)
    for i in range(ny):
        for j in range(nx):
            sl = np.s_[i * step:(i + 1) * step, j * step:(j + 1) * step]
            v = signal_image[sl][free[sl]]
            count[i, j] = v.size
            if v.size >= min_pixels:
                med[i, j] = np.median(v)
    return med, count


def cell_free_bkg_estimation(signal_image, free, step=ESTIMATION_STEP, smoothing_sigma=SMOOTHING_SIGMA,
                             min_pixels=MIN_PIXELS):
    """Background surface: tile medians of the cell-free pixels, smoothed with a Gaussian (sigma in px).
    Tiles with too few cell-free pixels are filled in by normalised smoothing of the others.
    -> (surface, tile medians, tile pixel counts)"""
    signal_image = np.asarray(signal_image, float)
    med, count = tile_medians(signal_image, free, step, min_pixels)
    valid = np.isfinite(med)
    if not valid.any():
        raise ValueError(f'no {step} px tile has {min_pixels} cell-free pixels: lower dilation, '
                         'estimation_step or min_pixels')
    H, W = signal_image.shape
    block = np.ones((step, step))
    fill = np.median(med[valid])
    if smoothing_sigma and smoothing_sigma > 0:
        values = np.kron(np.where(valid, med, 0.0), block)[:H, :W]
        weights = np.kron(valid.astype(float), block)[:H, :W]
        num = ndimage.gaussian_filter(values, smoothing_sigma)
        den = ndimage.gaussian_filter(weights, smoothing_sigma)
        surface = np.where(den > 1e-6, num / np.maximum(den, 1e-12), fill)
    else:
        surface = np.kron(np.where(valid, med, fill), block)[:H, :W]
    return surface, med, count


def noise_sd(signal_image, free):
    """Noise SD of the cell-free pixels from differences of horizontal and vertical neighbours:
    1.4826 x MAD / sqrt 2 (blind to slow structure such as glow or uneven background)."""
    s = np.asarray(signal_image, float)
    d = np.concatenate([(s[:, 1:] - s[:, :-1])[free[:, 1:] & free[:, :-1]],
                        (s[1:, :] - s[:-1, :])[free[1:, :] & free[:-1, :]]])
    if d.size < 10:
        return np.nan
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


# ---------------------------------------------------------------- subtraction
def back_sub(signal_image, phase_mask=None, dilation=DILATION, estimation_step=ESTIMATION_STEP,
             smoothing_sigma=SMOOTHING_SIGMA, show=False, cell_masks=None, camera_offset=0.0,
             flat_norm='median', min_pixels=MIN_PIXELS, return_info=False):
    """
    Subtracts the cell-free background surface from a fluorescence image.

    Parameters
    ----------
    signal_image: 2D numpy array - the image to be corrected
    phase_mask: boolean 2D array - cells, e.g. from get_inverted_mask(phase_image); optional if cell_masks is given
    dilation: px (Euclidean) around the cells excluded from the background (default 30)
    estimation_step: px, side of the square tiles in which the cell-free median is taken (default 128; it does
        not have to divide the image size)
    smoothing_sigma: px, Gaussian smoothing of the tile medians (default 100)
    show: True to plot the excluded pixels, the surface and the corrected image
    cell_masks: 2D array (labels or boolean) of more cells to exclude, e.g. the Omnipose masks
    camera_offset: counts the camera adds without light, removed from the flat-field profile (default 0)
    flat_norm: 'median' (default) or 'max' (version 1): normalisation of the flat-field profile
    min_pixels: cell-free pixels a tile needs to count (default 50)
    return_info: True to also return a dict with the noise SD and diagnostics

    Returns
    -------
    [0] the background-subtracted image divided by the flat-field profile (uneven excitation corrected)
    [1] the background-subtracted image
    [2] the background surface
    [3] the excluded pixels (True = cell or within `dilation` px of a cell)
    [4] (return_info) dict: noise_sd, background_median, surface_min / max, flat_min / max, free_fraction,
        n_tiles, n_empty_tiles, dilation
    """
    signal = np.asarray(signal_image, float)
    if phase_mask is None and cell_masks is None:
        raise ValueError('give phase_mask and / or cell_masks')
    cells = np.zeros(signal.shape, bool)
    if phase_mask is not None:
        cells |= np.asarray(phase_mask, bool)
    if cell_masks is not None:
        cells |= np.asarray(cell_masks) > 0
    free = cell_free_pixels(cells, dilation)
    surface, med, count = cell_free_bkg_estimation(signal, free, estimation_step, smoothing_sigma, min_pixels)
    bkg_cor_2 = signal - surface
    profile = surface - camera_offset
    if np.min(profile) <= 0:
        raise ValueError(f'camera_offset ({camera_offset}) is not below the background everywhere '
                         f'(surface minimum {surface.min():.1f})')
    if flat_norm == 'median':
        ref = np.median(profile)
    elif flat_norm == 'max':
        ref = np.max(profile)
    else:
        raise ValueError("flat_norm must be 'median' or 'max'")
    flat = profile / ref
    bkg_cor = bkg_cor_2 / flat
    excluded = ~free

    if show:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(18, 6), layout='constrained')
        v = signal[free]
        lo, hi = np.percentile(v, [1, 99]) if v.size else (signal.min(), signal.max())
        ax[0].imshow(np.where(free, signal, np.nan), cmap='magma', vmin=lo, vmax=hi)
        ax[0].set_title(f'cell-free pixels (> {dilation} px from any cell)')
        im = ax[1].imshow(surface, cmap='Blues')
        fig.colorbar(im, ax=ax[1], shrink=0.7)
        ax[1].set_title('background surface')
        ax[2].imshow(bkg_cor, cmap='Blues', vmin=0, vmax=np.percentile(bkg_cor, 99.9))
        ax[2].set_title('corrected image [0]')
        for a in ax:
            a.axis('off')
        plt.show()

    out = (bkg_cor, bkg_cor_2, surface, excluded)
    if not return_info:
        return out
    info = dict(noise_sd=noise_sd(signal, free), background_median=float(np.median(surface)),
                surface_min=float(surface.min()), surface_max=float(surface.max()),
                flat_min=float(flat.min()), flat_max=float(flat.max()), free_fraction=float(free.mean()),
                n_tiles=int(med.size), n_empty_tiles=int((~np.isfinite(med)).sum()), dilation=dilation)
    return out + (info,)


def background_pixels(corrected_image, excluded):
    """The cell-free pixels of a corrected image (e.g. for the residual background mean / SD)."""
    return np.asarray(corrected_image)[~np.asarray(excluded, bool)]
