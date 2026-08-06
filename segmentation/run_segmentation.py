#!/usr/bin/env python
"""
run_segmentation.py
====================
CLI entry point for the PhenoSuite segmentation module. Routes a single
multi-channel image through one of four cell-segmentation backends:

  sam       Segment Anything (Meta) — promptless automatic mask generation
  mesmer    DeepCell Mesmer         — nuclear + whole-cell segmentation
  stardist  StarDist                — star-convex nuclei segmentation
  cellpose  Cellpose                — generalist nuclei / cytoplasm segmentation (default)

Usage:
  python run_segmentation.py --image=PATH --out-dir=DIR [options]

Run with --help for the full option reference. Fully standalone-runnable —
no SLURM required for a single image.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import segmentation_utils as su

EPILOG = """
Sentinels that all mean "not set" for --label, --nuclear-channel,
--membrane-channel, --sam-checkpoint, and per-sample batch-inputs rows:
NULL, null, NA, none, None, '' (empty). Note '0' is NOT a sentinel here —
unlike other phenosuite-CLI modules, it's a legitimate 0-based channel
index. That makes every row of every batch-inputs/*.txt file safe to pad
with NULL (never with 0).

Outputs (written to --out-dir, prefixed with --label):
  {label}_mask.tif          uint32 label mask, same H x W as the input image
  {label}_centroids.csv     label, y, x, area — one row per segmented object.
                            Having label/y/x columns present matches the bare
                            "Mesmer" segmentation fingerprint RunPhenomenalist
                            auto-detects, so this file (plus any marker
                            intensities you join onto it) drops straight into
                            RunPhenomenalist / masquerade regardless of which
                            route produced it.
  {label}_overlay.png       nuclear channel with mask boundaries in red (QC)
  {label}_provenance.json   inputs, route, parameters, timing, git SHA

Route-specific notes:
  sam       Requires --sam-checkpoint (a downloaded ViT checkpoint .pth).
            No tiling and no cell-biology prior -- runs once on a pseudo-RGB
            composite of the nuclear (+ membrane) channel, so it doesn't
            scale to whole-slide QPTIFF/OME-TIFF images. Guarded by
            --sam-max-side (default 4096px on the longest side): larger
            images are refused with an error rather than silently running
            slow and low-quality. Crop to an ROI first, pick another route,
            or pass --sam-max-side=0 / a larger value to override.
  mesmer    Whole-image inference only; --tile-size is ignored.
            Needs --resolution (microns/pixel) for correct cell sizing.
  stardist  Nuclei only. Set --tile-size > 0 to use predict_instances_big
            for images too large to fit in memory whole.
  cellpose  Supports --membrane-channel for two-channel (cyto) segmentation;
            nuclear-only when omitted. Always tiles internally.

Environment variables:
  PHENOSUITE_GIT_SHA        recorded in provenance.json (else "unknown")
  PHENOSUITE_IMAGE_DIGEST   recorded in provenance.json (else "unknown")

Prerequisites:
  conda env create -f environment.yml && conda activate segmentation
  (trim environment.yml to the one or two routes you actually use — the
  four backends pull in both torch and tensorflow and are heavy together)
"""


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    required = parser.add_argument_group("required")
    required.add_argument("--image", required=True, help="multi-channel TIFF / OME-TIFF / QPTIFF path")
    required.add_argument("--out-dir", required=True, dest="out_dir", help="output directory (created if missing)")

    route = parser.add_argument_group("route")
    route.add_argument("--method", default="cellpose", choices=sorted(su.ROUTES),
                        help="segmentation backend to use [default: cellpose]")

    chan = parser.add_argument_group("channels")
    chan.add_argument("--nuclear-channel", dest="nuclear_channel", default=None,
                       help="channel name or 0-based index (e.g. DAPI) [default: auto-detect]")
    chan.add_argument("--membrane-channel", dest="membrane_channel", default=None,
                       help="channel name or 0-based index for whole-cell routes [default: nuclear-only]")

    params = parser.add_argument_group("route parameters")
    params.add_argument("--diameter", type=float, default=0.0,
                         help="expected object diameter in px for cellpose [default: 0 -> auto-estimate]")
    params.add_argument("--cellpose-model", dest="cellpose_model", default="cyto3",
                         help="cellpose model_type [default: cyto3]")
    params.add_argument("--stardist-model", dest="stardist_model", default="2D_versatile_fluo",
                         help="pretrained StarDist2D model name [default: 2D_versatile_fluo]")
    params.add_argument("--sam-checkpoint", dest="sam_checkpoint", default=None,
                         help="path to a SAM ViT checkpoint .pth (required for --method=sam)")
    params.add_argument("--sam-model-type", dest="sam_model_type", default="vit_b",
                         choices=["vit_h", "vit_l", "vit_b"],
                         help="SAM backbone matching --sam-checkpoint [default: vit_b]")
    params.add_argument("--sam-max-side", dest="sam_max_side", type=int, default=su.SAM_DEFAULT_MAX_SIDE,
                         help="safety guard: refuse --method=sam above this many px on the "
                              "longest side (SAM has no tiling / cell-biology prior, so it "
                              "doesn't scale to whole-slide images); 0 disables the guard "
                              f"[default: {su.SAM_DEFAULT_MAX_SIDE}]")
    params.add_argument("--resolution", type=float, default=0.5,
                         help="microns per pixel, used by mesmer for cell sizing [default: 0.5]")

    post = parser.add_argument_group("post-processing")
    post.add_argument("--min-size", dest="min_size", type=int, default=0,
                       help="drop objects smaller than this many pixels [default: 0 -> no filter]")
    post.add_argument("--tile-size", dest="tile_size", type=int, default=0,
                       help="tile size in px for stardist / cellpose on large images [default: 0 -> auto/whole-image]")
    post.add_argument("--tile-overlap", dest="tile_overlap", type=int, default=64,
                       help="tile overlap in px (stardist only) [default: 64]")

    misc = parser.add_argument_group("misc")
    misc.add_argument("--label", default=None, help="output filename prefix [default: image basename]")
    gpu_group = misc.add_mutually_exclusive_group()
    gpu_group.add_argument("--gpu", dest="gpu", action="store_true", default=None, help="force GPU")
    gpu_group.add_argument("--no-gpu", dest="gpu", action="store_false", help="force CPU")
    misc.add_argument("--no-overlay", dest="export_overlay", action="store_false", default=True,
                       help="skip the QC overlay PNG")
    misc.add_argument("--no-centroids", dest="export_centroids", action="store_false", default=True,
                       help="skip the centroids CSV")
    misc.add_argument("--seed", type=int, default=42, help="RNG seed [default: 42]")

    return parser


def clean_sentinel(value):
    return None if su.is_unset(value) else value


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    args.nuclear_channel = clean_sentinel(args.nuclear_channel)
    args.membrane_channel = clean_sentinel(args.membrane_channel)
    args.label = clean_sentinel(args.label)
    args.sam_checkpoint = clean_sentinel(args.sam_checkpoint)

    if not os.path.exists(args.image):
        parser.error(f"--image not found: {args.image}")

    if args.label is None:
        args.label = os.path.splitext(os.path.basename(args.image))[0]
        args.label = args.label[:-4] if args.label.lower().endswith(".ome") else args.label

    if args.method == "sam" and not args.sam_checkpoint:
        parser.error("--sam-checkpoint is required when --method=sam")

    su.run_segmentation(
        image_path=args.image,
        out_dir=args.out_dir,
        label=args.label,
        method=args.method,
        nuclear_channel=args.nuclear_channel,
        membrane_channel=args.membrane_channel,
        diameter=args.diameter,
        cellpose_model=args.cellpose_model,
        stardist_model=args.stardist_model,
        sam_checkpoint=args.sam_checkpoint,
        sam_model_type=args.sam_model_type,
        sam_max_side=args.sam_max_side,
        resolution=args.resolution,
        min_size=args.min_size,
        tile_size=args.tile_size,
        tile_overlap=args.tile_overlap,
        gpu=args.gpu,
        export_overlay=args.export_overlay,
        export_centroids=args.export_centroids,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
