#!/usr/bin/env python
"""
segmentation_utils.py
======================
Core image I/O, channel resolution, and per-route segmentation runners for
the PhenoSuite segmentation CLI.

Four interchangeable routes share one input/output contract:

  sam       Segment Anything (Meta) — promptless automatic mask generation
  mesmer    DeepCell Mesmer         — nuclear + whole-cell segmentation
  stardist  StarDist                — star-convex nuclei segmentation
  cellpose  Cellpose                — generalist nuclei / cytoplasm segmentation

Sourced by run_segmentation.py. Only numpy / pandas / tifffile / scikit-image /
matplotlib are imported at module load time — the heavy, route-specific
backends (torch, tensorflow, deepcell, stardist, cellpose, segment-anything)
are imported lazily inside their respective run_* function so that `--help`
and unrelated routes work without every backend installed.
"""

from __future__ import annotations

import json
import os
import sys
import time
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import tifffile
from tifffile import TiffFile

SENTINELS = {"", "null", "none", "na"}


def is_unset(value) -> bool:
    """True for the batch-inputs NULL-padding sentinels (NULL/null/NA/none/
    None/''). Deliberately does NOT include '0' — unlike the numeric-only
    sentinel convention used elsewhere in this repo, '0' is a legitimate
    0-based channel index here (often the DAPI/nuclear channel)."""
    if value is None:
        return True
    return str(value).strip().lower() in SENTINELS


# ---------------------------------------------------------------------------
# Image I/O — channel-aware TIFF / OME-TIFF / QPTIFF loading
# ---------------------------------------------------------------------------

def detect_tiff_type(path: str) -> str:
    with TiffFile(path) as tif:
        if tif.is_ome:
            return "ome"
        try:
            first_page = tif.series[0].pages[0]
            desc = getattr(first_page, "description", "") or ""
            if "Biomarker" in desc:
                return "qptiff"
        except Exception:
            pass
        return "standard"


def _channel_names_ome(tif: TiffFile) -> list:
    try:
        root = ElementTree.fromstring(tif.ome_metadata)
        for ns_uri in (
            "http://www.openmicroscopy.org/Schemas/OME/2016-06",
            "http://www.openmicroscopy.org/Schemas/OME/2015-01",
            "http://www.openmicroscopy.org/Schemas/OME/2013-06",
        ):
            channels = root.findall(".//ome:Channel", {"ome": ns_uri})
            if channels:
                names = [ch.get("Name") for ch in channels]
                if all(names):
                    return names
    except Exception:
        pass
    return [f"Channel_{i}" for i in range(len(tif.series[0].pages))]


def _channel_names_qptiff(tif: TiffFile) -> list:
    names = []
    for page in tif.series[0].pages:
        desc = getattr(page, "description", "") or ""
        try:
            el = ElementTree.fromstring(desc).find("Biomarker")
            names.append(el.text.replace(" ", "-") if el is not None else f"Channel_{page.index}")
        except Exception:
            names.append(f"Channel_{page.index}")
    return names


def get_channel_names(path: str) -> list:
    """Return channel names by reading only TIFF/OME-XML/QPTIFF metadata —
    no pixel data. Auto-detects OME-TIFF / QPTIFF / plain-TIFF channel
    naming, reusing the masquerade detection conventions so channel names
    line up across modules. Multiplex TIFFs store one channel per page
    (SamplesPerPixel=1, non-interleaved), so this never decodes a full,
    possibly dozens-of-channels-wide image just to list its channel names."""
    tiff_type = detect_tiff_type(path)
    with TiffFile(path) as tif:
        n_pages = len(tif.series[0].pages) if tif.series[0].pages else 1
        if tiff_type == "ome":
            names = _channel_names_ome(tif)
        elif tiff_type == "qptiff":
            names = _channel_names_qptiff(tif)
        else:
            names = [f"Channel_{i}" for i in range(n_pages)]

    if len(names) != n_pages:
        names = [f"Channel_{i}" for i in range(n_pages)]
    return names


def read_channel(path: str, index: int):
    """Decode a single channel page as float32, without touching the other
    channels. This is the memory-critical path for QPTIFF/OME-TIFF inputs
    that routinely carry 20-60+ marker channels but only need 1-2 decoded
    (nuclear ± membrane) for segmentation."""
    with TiffFile(path) as tif:
        pages = tif.series[0].pages
        if index < 0 or index >= len(pages):
            raise ValueError(f"channel index {index} out of range (0..{len(pages) - 1})")
        arr = pages[index].asarray()

    if arr.ndim > 2:
        # Rare: a single page unexpectedly packs multiple samples/channels
        # (e.g. interleaved RGB) instead of the usual one-page-per-channel
        # layout. Fall back to the first band rather than decoding the rest
        # of the file.
        arr = arr[..., 0] if arr.shape[-1] <= 8 else arr[0]

    return arr.astype(np.float32)


def resolve_channel_index(channel_names, spec, fallback_patterns=()):
    """Resolve a channel spec (0-based index, exact/partial name, or unset)
    to an (index, resolved name) pair, using only the channel-name list —
    no pixel data required. fallback_patterns are regexes tried in order
    when spec is unset, e.g. ('dapi', 'hoechst') for nuclear auto-detection."""
    if is_unset(spec):
        import re
        for pat in fallback_patterns:
            for i, name in enumerate(channel_names):
                if re.search(pat, name, re.IGNORECASE):
                    return i, name
        return None, None

    spec = str(spec).strip()
    if spec.lstrip("-").isdigit():
        idx = int(spec)
        if idx < 0 or idx >= len(channel_names):
            raise ValueError(f"channel index {idx} out of range (0..{len(channel_names) - 1})")
        return idx, channel_names[idx]

    for i, name in enumerate(channel_names):
        if name.lower() == spec.lower():
            return i, name
    for i, name in enumerate(channel_names):
        if spec.lower() in name.lower():
            return i, name
    raise ValueError(f"channel '{spec}' not found among {channel_names}")


def normalize_percentile(img, low=1.0, high=99.8):
    lo, hi = np.percentile(img, (low, high))
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((img - lo) / (hi - lo), 0, 1).astype(np.float32)


def make_pseudo_rgb(nuclear, membrane=None):
    """Build a uint8 3-channel image for models (SAM) that expect RGB input."""
    n = (normalize_percentile(nuclear) * 255).astype(np.uint8)
    m = (normalize_percentile(membrane) * 255).astype(np.uint8) if membrane is not None else np.zeros_like(n)
    z = np.zeros_like(n)
    return np.stack([m, n, z], axis=-1)


def resolve_gpu(requested):
    """requested: True/False/None. None means auto-detect via torch (falls
    back to CPU silently if torch or a GPU isn't available)."""
    if requested is not None:
        return bool(requested)
    try:
        import torch
        return bool(torch.cuda.is_available())
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Route runners — each returns a uint32 label mask, shape (H, W)
# ---------------------------------------------------------------------------

def run_cellpose(nuclear, membrane, model_name, diameter, gpu, tile_size):
    try:
        from cellpose import models
    except ImportError as exc:
        raise RuntimeError(
            "--method=cellpose requires the 'cellpose' package.\n"
            "  Install it with:  pip install cellpose\n"
            "  (see segmentation/environment.yml)"
        ) from exc

    model = models.Cellpose(gpu=gpu, model_type=model_name)
    if membrane is not None:
        img = np.stack([membrane, nuclear], axis=-1)
        channels = [1, 2]
    else:
        img = nuclear
        channels = [0, 0]

    eval_kwargs = dict(diameter=(diameter if diameter and diameter > 0 else None),
                        channels=channels, tile=True)
    if tile_size and tile_size > 0:
        eval_kwargs["bsize"] = tile_size

    masks, _flows, _styles, _diams = model.eval(img, **eval_kwargs)
    return masks.astype(np.uint32)


def run_stardist(nuclear, model_name, tile_size, tile_overlap):
    try:
        from csbdeep.utils import normalize as _cs_normalize
        from stardist.models import StarDist2D
    except ImportError as exc:
        raise RuntimeError(
            "--method=stardist requires the 'stardist' package (+ tensorflow).\n"
            "  Install it with:  pip install stardist tensorflow\n"
            "  (see segmentation/environment.yml)"
        ) from exc

    model = StarDist2D.from_pretrained(model_name)
    img = _cs_normalize(nuclear, 1, 99.8)

    if tile_size and tile_size > 0:
        labels, _details = model.predict_instances_big(
            img, axes="YX", block_size=tile_size,
            min_overlap=tile_overlap, context=tile_overlap,
        )
    else:
        labels, _details = model.predict_instances(img)
    return labels.astype(np.uint32)


def run_mesmer(nuclear, membrane, resolution):
    try:
        from deepcell.applications import Mesmer
    except ImportError as exc:
        raise RuntimeError(
            "--method=mesmer requires the 'deepcell' package (+ tensorflow).\n"
            "  Install it with:  pip install deepcell tensorflow\n"
            "  (see segmentation/environment.yml)"
        ) from exc

    membrane_ch = membrane if membrane is not None else nuclear
    compartment = "whole-cell" if membrane is not None else "nuclear"
    stack = np.stack([nuclear, membrane_ch], axis=-1)[np.newaxis, ...]
    app = Mesmer()
    labeled = app.predict(stack, image_mpp=resolution, compartment=compartment)
    return labeled[0, ..., 0].astype(np.uint32)


SAM_DEFAULT_MAX_SIDE = 4096


def run_sam(nuclear, membrane, checkpoint, model_type, gpu, max_side=SAM_DEFAULT_MAX_SIDE):
    """max_side is a safety guard, not a quality setting: SamAutomaticMaskGenerator
    has no tiling and no cell-biology prior, so running it on a whole-slide
    QPTIFF/OME-TIFF (tens of thousands of px per side) is both slow and
    low-quality compared to cellpose/mesmer/stardist at that scale. Refuses
    above max_side px on the longest side unless max_side<=0 (no limit) --
    crop to an ROI first, or pass a larger/0 --sam-max-side to override."""
    longest_side = max(nuclear.shape)
    if max_side and max_side > 0 and longest_side > max_side:
        raise ValueError(
            f"--method=sam: image is {nuclear.shape[0]}x{nuclear.shape[1]}px "
            f"(longest side {longest_side} > --sam-max-side={max_side}). SAM has no "
            "tiling and no cell-biology prior, so it doesn't scale to whole-slide "
            "QPTIFF/OME-TIFF images -- crop to a region of interest first, or use "
            "--method=cellpose / mesmer / stardist for full images. To override "
            "this guard anyway, pass --sam-max-side=0 (no limit) or a larger value."
        )

    try:
        import torch
        from segment_anything import SamAutomaticMaskGenerator, sam_model_registry
    except ImportError as exc:
        raise RuntimeError(
            "--method=sam requires 'torch' and 'segment-anything'.\n"
            "  Install with:  pip install torch torchvision "
            "git+https://github.com/facebookresearch/segment-anything.git\n"
            "  (see segmentation/environment.yml)"
        ) from exc

    if not checkpoint or not os.path.exists(checkpoint):
        raise FileNotFoundError(
            f"--sam-checkpoint not found: '{checkpoint}'. Download a SAM ViT "
            "checkpoint (vit_h/vit_l/vit_b) and pass its path explicitly."
        )

    device = "cuda" if gpu and torch.cuda.is_available() else "cpu"
    sam = sam_model_registry[model_type](checkpoint=checkpoint)
    sam.to(device)
    generator = SamAutomaticMaskGenerator(sam)

    rgb = make_pseudo_rgb(nuclear, membrane)
    results = generator.generate(rgb)

    label_mask = np.zeros(rgb.shape[:2], dtype=np.uint32)
    for i, r in enumerate(sorted(results, key=lambda r: r["area"], reverse=True), start=1):
        label_mask[r["segmentation"]] = i
    return label_mask


ROUTES = {"sam", "mesmer", "stardist", "cellpose"}


# ---------------------------------------------------------------------------
# Post-processing + export
# ---------------------------------------------------------------------------

def filter_min_size(label_mask, min_size):
    if not min_size or min_size <= 0:
        return label_mask
    from skimage.morphology import remove_small_objects
    return remove_small_objects(label_mask, min_size=min_size).astype(np.uint32)


def extract_centroids(label_mask):
    """Regionprops -> a `label, y, x, area` table. Having `label`, `y`, `x`
    columns present (order doesn't matter) matches the bare 'Mesmer'
    segmentation-table fingerprint RunPhenomenalist auto-detects
    (phenomenalist-utils.R: mesmer_base <- c("label","y","x")), so the
    centroid CSV from any of the four routes drops straight into
    RunPhenomenalist / masquerade without extra reformatting."""
    from skimage.measure import regionprops_table

    props = regionprops_table(label_mask, properties=("label", "centroid", "area"))
    df = pd.DataFrame({
        "label": props["label"],
        "y": props["centroid-0"],
        "x": props["centroid-1"],
        "area": props["area"],
    })
    return df[["label", "y", "x", "area"]]


def write_overlay(nuclear, label_mask, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    boundaries = find_boundaries(label_mask, mode="inner")
    base = normalize_percentile(nuclear)
    rgb = np.stack([base, base, base], axis=-1)
    rgb[boundaries] = [1.0, 0.0, 0.0]
    plt.imsave(out_path, rgb)


def write_provenance(path, **fields):
    payload = {
        "tool": "phenosuite-segmentation-CLI",
        "run_date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python_version": sys.version.split()[0],
        "git_sha": os.environ.get("PHENOSUITE_GIT_SHA", "unknown"),
        "image_digest": os.environ.get("PHENOSUITE_IMAGE_DIGEST", "unknown"),
    }
    payload.update(fields)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)


def run_segmentation(
    image_path,
    out_dir,
    label,
    method,
    nuclear_channel=None,
    membrane_channel=None,
    diameter=0.0,
    cellpose_model="cyto3",
    stardist_model="2D_versatile_fluo",
    sam_checkpoint=None,
    sam_model_type="vit_b",
    sam_max_side=SAM_DEFAULT_MAX_SIDE,
    resolution=0.5,
    min_size=0,
    tile_size=0,
    tile_overlap=64,
    gpu=None,
    export_overlay=True,
    export_centroids=True,
    seed=42,
):
    if method not in ROUTES:
        raise ValueError(f"--method must be one of {sorted(ROUTES)} (got '{method}')")

    np.random.seed(seed)
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    # Only decode the 1-2 channels actually needed for segmentation, not the
    # full stack -- QPTIFF/OME-TIFF multiplex panels routinely carry 20-60+
    # marker channels, and decoding all of them just to segment on nuclear
    # (+ membrane) wastes an order of magnitude of memory on large images.
    channel_names = get_channel_names(image_path)

    nuclear_idx, nuclear_name = resolve_channel_index(
        channel_names, nuclear_channel,
        fallback_patterns=(r"dapi", r"hoechst", r"nucle", r"nuclear"),
    )
    if nuclear_idx is None:
        raise ValueError(
            f"could not resolve a nuclear channel from {channel_names}; "
            "pass --nuclear-channel explicitly"
        )
    nuclear = read_channel(image_path, nuclear_idx)

    membrane_idx, membrane_name = resolve_channel_index(
        channel_names, membrane_channel,
        fallback_patterns=(r"membrane", r"cd45", r"na.?k.?atpase", r"e.?cadherin"),
    )
    membrane = read_channel(image_path, membrane_idx) if membrane_idx is not None else None

    gpu_resolved = resolve_gpu(gpu)

    print("=" * 60)
    print("SEGMENTATION CONFIGURATION")
    print("=" * 60)
    print(f"image            --> {image_path}")
    print(f"out dir          --> {out_dir}")
    print(f"label            --> {label}")
    print(f"method           --> {method}")
    print(f"channels         --> {channel_names}")
    print(f"nuclear channel  --> {nuclear_name}")
    print(f"membrane channel --> {membrane_name}")
    print(f"gpu              --> {gpu_resolved}")
    print("=" * 60)

    if method == "cellpose":
        mask = run_cellpose(nuclear, membrane, cellpose_model, diameter, gpu_resolved, tile_size)
    elif method == "stardist":
        mask = run_stardist(nuclear, stardist_model, tile_size, tile_overlap)
    elif method == "mesmer":
        if tile_size and tile_size > 0:
            print("note: --tile-size is ignored for method=mesmer (whole-image inference only)")
        mask = run_mesmer(nuclear, membrane, resolution)
    elif method == "sam":
        if tile_size and tile_size > 0:
            print("note: --tile-size is ignored for method=sam (whole-image inference only)")
        mask = run_sam(nuclear, membrane, sam_checkpoint, sam_model_type, gpu_resolved, sam_max_side)
    else:  # pragma: no cover - guarded above
        raise AssertionError(method)

    mask = filter_min_size(mask, min_size)
    n_cells = int(len(np.unique(mask)) - (1 if 0 in mask else 0))
    elapsed = time.time() - t0
    print(f"segmented {n_cells} objects in {elapsed:.1f}s")

    outputs = {}

    mask_path = os.path.join(out_dir, f"{label}_mask.tif")
    tifffile.imwrite(mask_path, mask, compression="zlib")
    outputs["mask"] = mask_path

    if export_centroids:
        centroids_path = os.path.join(out_dir, f"{label}_centroids.csv")
        extract_centroids(mask).to_csv(centroids_path, index=False)
        outputs["centroids"] = centroids_path

    if export_overlay:
        overlay_path = os.path.join(out_dir, f"{label}_overlay.png")
        write_overlay(nuclear, mask, overlay_path)
        outputs["overlay"] = overlay_path

    prov_path = os.path.join(out_dir, f"{label}_provenance.json")
    write_provenance(
        prov_path,
        image=os.path.abspath(image_path),
        label=label,
        method=method,
        n_cells=n_cells,
        elapsed_sec=round(elapsed, 2),
        channels=channel_names,
        nuclear_channel=nuclear_name,
        membrane_channel=membrane_name,
        diameter=diameter,
        cellpose_model=cellpose_model if method == "cellpose" else None,
        stardist_model=stardist_model if method == "stardist" else None,
        sam_model_type=sam_model_type if method == "sam" else None,
        sam_max_side=sam_max_side if method == "sam" else None,
        resolution=resolution if method == "mesmer" else None,
        min_size=min_size,
        tile_size=tile_size,
        tile_overlap=tile_overlap,
        gpu=gpu_resolved,
        seed=seed,
        outputs=outputs,
    )
    outputs["provenance"] = prov_path

    print(f"done. outputs written to: {out_dir}")
    return outputs
