#!/usr/bin/env python
"""
Masquerade Interface Script
Handles OME-TIFF and QPTIFF input formats.
Outputs .ome.tiff so QuPath uses the Bio-Formats OME reader,
which respects the manually constructed OME-XML with SamplesPerPixel=1.
"""

import sys
import re
import numpy as np

# Command line arguments:
# ${image} ${spatial_metadata} ${outPath} ${relevant_markers} ${adjust_coords} ${compress} ${radius} ${filled} ${num_points} ${preFilter_masks} ${target_size}
image_source    = sys.argv[1]
metadata        = sys.argv[2]
out_path        = sys.argv[3]
marker_metadata = sys.argv[4]
adjust_coords   = sys.argv[5]
compress        = sys.argv[6]
radius          = sys.argv[7]
filled          = sys.argv[8]
num_points      = sys.argv[9]
preFilter_masks = sys.argv[10]
target_size     = sys.argv[11]

from Masquerade import Masquerade

masquerade = Masquerade()

adjust_coords   = adjust_coords == 'True'
compress        = compress == 'True'
preFilter_masks = preFilter_masks == 'True'
filled          = filled == 'True'

if marker_metadata == 'None':
	marker_metadata = None

print('='*60)
print('MASQUERADE CONFIGURATION')
print('='*60)
print('image source     --> ' + image_source)
print('spatial annos    --> ' + metadata)
print('out path         --> ' + out_path)
print('marker metadata  --> ' + str(marker_metadata))
print('adjust coords:   ' + str(adjust_coords))
print('compress:        ' + str(compress))
print('radius:          ' + str(radius))
print('fill circle:     ' + str(filled))
print('num points:      ' + str(num_points))
print('pre-filter masks:' + str(preFilter_masks))
print('target size      --> ' + str(target_size))
print('='*60)

PYRAMIDAL_TILE_SIZE  = 512
PYRAMIDAL_COMPRESSION = None

masquerade.relevant_markers = marker_metadata
masquerade.spatial_anno     = metadata
masquerade.image_source     = image_source
masquerade.compress         = compress
masquerade.target_size      = float(target_size)
masquerade.filter_img       = preFilter_masks
masquerade.radius           = int(radius)
masquerade.filled           = filled
masquerade.num_points       = int(num_points)
masquerade.adjust_coords    = adjust_coords

print('\n' + '='*60)
print('PROCESSING STEPS')
print('='*60)

# Step 1: Detect TIFF type
print('\n[1/5] Detecting TIFF format...')
tiff_type = masquerade.detect_tiff_type()
print(f'      Format detected: {tiff_type.upper()}')

# Step 2: Pre-process image
print('\n[2/5] Pre-processing image...')
print(f'      Raw size before processing: {masquerade.raw_size} GB')

image, set_subset_x, set_subset_y, spatial_metadata = masquerade.PreProcessImage()

print(f'      Raw size after processing: {masquerade.raw_size} GB')
print(f'      Image shape: {image.shape}')
print(f'      Compression factor: {masquerade.compression_factor}')

# Step 3: Generate circle masks
print('\n[3/5] Generating circle masks...')

print(f'      Validating spatial annotations...')
print(f'      File: {metadata}')
import pandas as pd
spatial_check = pd.read_csv(metadata)
print(f'      Total rows: {len(spatial_check):,}')
print(f'      Columns: {spatial_check.columns.tolist()}')

if 'cluster' in spatial_check.columns:
	n_clusters    = spatial_check['cluster'].nunique()
	cluster_sizes = spatial_check.groupby('cluster').size()
	valid_clusters = (cluster_sizes > 1).sum()
	print(f'      Total unique clusters: {n_clusters}')
	print(f'      Clusters with >1 cell: {valid_clusters}')
	print(f'      Cluster size range: {cluster_sizes.min()} to {cluster_sizes.max()}')
	if valid_clusters == 0:
		print('      ⚠️  ERROR: No valid clusters (all have ≤1 cell)!')
else:
	print('      ⚠️  ERROR: No "cluster" column found in spatial annotations!')

channels = masquerade.get_circle_masks(
	image=image,
	metadata=spatial_metadata,
	set_subset_x=set_subset_x,
	set_subset_y=set_subset_y
)
print(f'      Generated {len(channels)} mask channels')

# Step 4: Process marker channels and prepare output path
print('\n[4/5] Processing marker channels...')

import os

def force_ome_tiff_ext(filename):
	"""Ensure filename ends with .ome.tiff regardless of input extension."""
	for ext in ('.ome.tiff', '.ome.tif', '.tiff', '.tif'):
		if filename.endswith(ext):
			return filename[:-len(ext)] + '.ome.tiff'
	return filename + '.ome.tiff'

# Handle output path
if os.path.isdir(out_path):
	out_dir = out_path
	image_base   = os.path.splitext(os.path.basename(image_source))[0]
	out_filename = f"{image_base}_masks.ome.tiff"
	print(f'      Output is a directory, auto-generating filename: {out_filename}')

elif out_path.endswith('/') or (not os.path.splitext(out_path)[1]):
	out_dir = out_path.rstrip('/')
	image_base   = os.path.splitext(os.path.basename(image_source))[0]
	out_filename = f"{image_base}_masks.ome.tiff"
	print(f'      Output path treated as directory, auto-generating filename: {out_filename}')

else:
	out_dir      = os.path.dirname(out_path)
	out_filename = force_ome_tiff_ext(os.path.basename(out_path))
	print(f'      Output is a file path: {out_filename}')

if compress:
	print(f'      Compressing masks to {masquerade.target_size} GB')
	print(f'      Using compression factor: {masquerade.compression_factor}')

	channels, channel_names = masquerade.compress_marker_channels(
		channels,
		spatial_metadata=spatial_metadata
	)

	print(f'      Total channels after compression: {len(channels)}')
	print(f'      Channel names: {channel_names}')
	print('      Channel shapes:')
	for k in channels.keys():
		print(f'        {k}: {channels[k].shape}')

	# Insert compression factor before .ome.tiff
	base = out_filename[:-len('.ome.tiff')] if out_filename.endswith('.ome.tiff') else out_filename
	out_filename = f"{base}_cf{np.round(masquerade.compression_factor, 3)}.ome.tiff"
	print(f'      Updated filename with compression factor: {out_filename}')

else:
	print('      No compression requested')
	channels, channel_names = masquerade.get_continuous_channels(
		channels,
		set_subset_x=set_subset_x,
		set_subset_y=set_subset_y
	)
	print(f'      Total channels: {len(channels)}')
	print(f'      Channel names: {channel_names}')

out_path = os.path.join(out_dir, out_filename)
print(f'      Final output path: {out_path}')

# Step 5: Write output
print('\n[5/5] Writing OME-BigTIFF output...')
print(f'      Output file: {out_path}')
print(f'      Number of channels: {len(channels)}')

out_dir_final = os.path.dirname(out_path)
if out_dir_final and not os.path.exists(out_dir_final):
	print(f'      Creating output directory: {out_dir_final}')
	os.makedirs(out_dir_final, exist_ok=True)

masquerade.write_ome_bigTiff(
	channels=channels,
	out=out_path,
	channels_to_keep=list(channels.keys()),
	tile_size=PYRAMIDAL_TILE_SIZE,
	compression=PYRAMIDAL_COMPRESSION,
)

print('\n' + '='*60)
print('PROCESSING COMPLETE')
print('='*60)
print(f'Output saved to: {out_path}')
print(f'Total channels: {len(channels)}')
print(f'Final size: ~{masquerade.target_size if compress else masquerade.raw_size} GB')
print('='*60 + '\n')
