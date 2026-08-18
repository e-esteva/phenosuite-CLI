import imagecodecs
import skimage.io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tifffile
import sys
import xml
from tifffile import imread, TiffFile
from xml.etree import ElementTree
import os
import re
from scipy import ndimage
from scipy.ndimage import zoom
import tqdm

class Masquerade:
	def __init__(self):
		self.image_source=''
		self.spatial_anno=''
		self.raw_size=0
		self.target_size=4
		self.adjust_coords=True
		self.filled=True
		self.radius=10
		self.num_points=100
		self.compress=False
		self.relevant_markers=None
		self.mapping=None
		self.skip_classes=None
		self.filter_img=False
		self.compression_factor=1
		self.tiff_type=None  # Will be detected: 'ome', 'qptiff', or 'standard'

	def pullExt(self, file):
		return re.split('[.]',file)[-1]

	def detect_tiff_type(self):
		"""
		Detect the type of TIFF file (OME-TIFF, QPTIFF, or standard TIFF)
		"""
		with TiffFile(self.image_source) as tif:
			# Check if it's an OME-TIFF
			if tif.is_ome:
				print("Detected OME-TIFF format")
				return 'ome'

			# Check if it's a QPTIFF by looking for Biomarker tags in description
			try:
				first_page = tif.series[0].pages[0]
				if hasattr(first_page, 'description') and first_page.description:
					if 'Biomarker' in first_page.description:
						print("Detected QPTIFF format")
						return 'qptiff'
			except:
				pass

			print("Detected standard TIFF format")
			return 'standard'

	def get_channel_names_ome(self, tif):
		"""
		Extract channel names from OME-TIFF metadata
		"""
		channel_names = []

		try:
			ome_metadata = tif.ome_metadata
			if ome_metadata:
				root = ElementTree.fromstring(ome_metadata)
				ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

				# Try different namespace versions
				if not root.findall('.//ome:Channel', ns):
					ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2015-01'}
				if not root.findall('.//ome:Channel', ns):
					ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2013-06'}

				channels = root.findall('.//ome:Channel', ns)
				for ch in channels:
					name = ch.get('Name')
					if name:
						channel_names.append(name)

				if channel_names:
					return channel_names
		except Exception as e:
			print(f"Could not parse OME metadata: {e}")

		# Fallback
		if not channel_names:
			for i, page in enumerate(tif.series[0].pages):
				channel_names.append(f"Channel_{i}")

		return channel_names

	def get_channel_names_qptiff(self, tif):
		"""
		Extract channel names from QPTIFF metadata
		"""
		channel_names = []

		for page in tif.series[0].pages:
			if hasattr(page, 'description') and page.description:
				try:
					page_name_element = ElementTree.fromstring(page.description).find('Biomarker')
					if page_name_element is not None:
						page_name = page_name_element.text
						page_name = page_name.replace(' ', '-')
						channel_names.append(page_name)
					else:
						channel_names.append(f"Channel_{page.index}")
				except:
					channel_names.append(f"Channel_{page.index}")
			else:
				channel_names.append(f"Channel_{page.index}")

		return channel_names

	def PreProcessImage(self):
		# Detect TIFF type
		self.tiff_type = self.detect_tiff_type()

		# import spatial metadata:
		spatial_metadata=pd.read_csv(self.spatial_anno)

		# import qptiff:
		image = skimage.io.imread(fname=str(self.image_source),plugin='tifffile')

		subset_x=np.array(spatial_metadata['x']).astype(int)
		set_subset_x=np.array([x for x in set(subset_x)]).astype(int)

		subset_y=np.array(spatial_metadata['y']).astype(int)
		set_subset_y=np.array([x for x in set(subset_y)]).astype(int)

		if self.adjust_coords:
			image=image[:,np.min(set_subset_y):np.max(set_subset_y),:]
			image=image[:,:,np.min(set_subset_x):np.max(set_subset_x)]

		# compute raw image size as uint8 array in GB:
		raw_img_size=(np.array(image,dtype='uint8').nbytes/1e9)

		self.raw_size = raw_img_size

		return image,set_subset_x,set_subset_y,spatial_metadata

	def compress_marker_channels(self, channels, spatial_metadata):
		"""
		Updated to handle both OME-TIFF and QPTIFF formats
		"""
		if self.tiff_type is None:
			self.tiff_type = self.detect_tiff_type()

		channel_names=[]

		with TiffFile(self.image_source) as tif:
			if self.tiff_type == 'ome':
				all_channel_names = self.get_channel_names_ome(tif)
			elif self.tiff_type == 'qptiff':
				all_channel_names = self.get_channel_names_qptiff(tif)
			else:
				all_channel_names = [f"Channel_{i}" for i in range(len(tif.series[0].pages))]

			print(f"Found {len(all_channel_names)} channels")

			for page_idx, page in enumerate(tif.series[0].pages):
				page_name = all_channel_names[page_idx] if page_idx < len(all_channel_names) else f"Channel_{page_idx}"

				print(f"{page_idx}\t{page_name}")

				# Handle duplicate channel names
				if page_name in channels:
					page_name = page_name + '_'+str(channel_names.count(page_name)+1)
				channel_names.append(page_name)

				should_include = True
				if self.relevant_markers is not None:
					relevant_markers_ = pd.read_csv(self.relevant_markers)
					markers_use = relevant_markers_['x']
					markers_mod = ["".join(re.split('_',x)) for x in markers_use]
					markers_mod_tmp = [re.sub('_','-',x) for x in markers_use]
					markers_mod.extend(markers_mod_tmp)

					channel_mod = "".join(re.split("-",page_name))
					should_include = (page_name in markers_use or page_name in markers_mod)

				if should_include:
					page_tmp = page.asarray()

					if self.adjust_coords:
						subset_x = np.array(spatial_metadata['x']).astype(int)
						set_subset_x = np.array([x for x in set(subset_x)]).astype(int)

						subset_y = np.array(spatial_metadata['y']).astype(int)
						set_subset_y = np.array([x for x in set(subset_y)]).astype(int)

						page_tmp = page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
						page_tmp = page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]

					if self.compress:
						page_tmp = zoom(page_tmp, self.compression_factor)

					channels[page_name] = page_tmp

		return channels, channel_names

	@staticmethod
	def writeMaskTiff(channels, outPath):
		import os

		num_channels=len(channels.keys())
		print('Markers retained: ')
		print(channels.keys())

		out_dir = os.path.dirname(outPath)
		if out_dir and not os.path.exists(out_dir):
			print(f'Creating output directory: {out_dir}')
			os.makedirs(out_dir, exist_ok=True)

		new_tiff=np.stack([channels[x] for x in channels.keys()])
		new_tiff=np.array(new_tiff,dtype='uint8')

		print("writing new tiff")
		channels=list(channels.keys())
		ijmetadata={"Labels":channels}

		tifffile.imwrite(str(outPath), new_tiff, imagej=True, metadata=ijmetadata)

	@staticmethod
	def generate_circle(x_center, y_center, radius, num_points=100):
		"""
		Generate points for a filled circle.

		Parameters:
		-----------
		x_center : float
			x-coordinate of circle center
		y_center : float
			y-coordinate of circle center
		radius : float
			radius of the circle
		num_points : int, optional
			number of points to generate (default: 100)

		Returns:
		--------
		tuple
			Arrays of x and y coordinates of points in the circle
		"""
		theta = np.linspace(0, 2*np.pi, num_points)
		x = x_center + radius * np.cos(theta)
		y = y_center + radius * np.sin(theta)
		return x, y

	@staticmethod
	def generate_filled_circle(cx, cy, r):
		coordinates = []
		for x in range(cx - r, cx + r + 1):
			for y in range(cy - r, cy + r + 1):
				if (x - cx)**2 + (y - cy)**2 <= r**2:
					coordinates.append((x, y))
		return coordinates

	def get_circle_masks(self, image, metadata, set_subset_x, set_subset_y):
		binned_metadata=[metadata[metadata['cluster']==x] for x in list(set(metadata['cluster']))]

		print(f"      Total clusters to process: {len(binned_metadata)}")
		print(f"      Image shape for masks: {image.shape[1:3]}")

		masks = {}
		masks_preFiltered = {}
		channels = {}

		valid_clusters = 0
		skipped_clusters = 0

		for cluster in tqdm.tqdm(binned_metadata, desc="Processing clusters"):
			if cluster.shape[0] > 1:
				valid_clusters += 1
				cluster_id=cluster['cluster'][cluster['cluster'].keys()[0]]

				if self.adjust_coords:
					cluster_x=cluster['x']-np.min(set_subset_x)
					cluster_y=cluster['y']-np.min(set_subset_y)
				else:
					cluster_x=cluster['x']
					cluster_y=cluster['y']

				if valid_clusters <= 3 or valid_clusters % 10 == 0:
					print(f'\nProcessing cluster: {cluster_id} (#{valid_clusters}, {cluster.shape[0]} cells)')

				spatial_data_cluster_x=np.array([int(x) for x in cluster_x])
				spatial_data_cluster_y=np.array([int(x) for x in cluster_y])

				image_tmp=image
				cluster_channel=str(cluster_id)+'_mask-original'

				spatial_data_cluster_x0=[int(x) for x in cluster_x]
				spatial_data_cluster_y0=[int(x) for x in cluster_y]

				if self.filled:
					circle_coords = [self.generate_filled_circle(cx=spatial_data_cluster_x0[x], cy=spatial_data_cluster_y0[x], r=self.radius) for x in range(len(spatial_data_cluster_x0))]
					spatial_data_cluster = [y for x in circle_coords for y in x]
					spatial_data_cluster_x = [list(x)[0] for x in spatial_data_cluster]
					spatial_data_cluster_y = [list(x)[1] for x in spatial_data_cluster]

				else:
					circle_coords = [self.generate_circle(x_center=spatial_data_cluster_x0[x], y_center=spatial_data_cluster_y0[x], radius=self.radius, num_points=self.num_points) for x in range(len(spatial_data_cluster_x0))]
					spatial_data_cluster_x_ = [x[0] for x in circle_coords]
					spatial_data_cluster_x = [y for x in spatial_data_cluster_x_ for y in x]
					spatial_data_cluster_y_ = [x[1] for x in circle_coords]
					spatial_data_cluster_y = [y for x in spatial_data_cluster_y_ for y in x]

					spatial_data_cluster_x.extend(spatial_data_cluster_x0)
					spatial_data_cluster_x1=[int(x)+1 for x in spatial_data_cluster_x0]
					spatial_data_cluster_x2=[int(x)+2 for x in spatial_data_cluster_x0]
					spatial_data_cluster_x3=[int(x)-1 for x in spatial_data_cluster_x0]
					spatial_data_cluster_x4=[int(x)-2 for x in spatial_data_cluster_x0]
					spatial_data_cluster_x.extend(spatial_data_cluster_x1)
					spatial_data_cluster_x.extend(spatial_data_cluster_x2)
					spatial_data_cluster_x.extend(spatial_data_cluster_x3)
					spatial_data_cluster_x.extend(spatial_data_cluster_x4)

					spatial_data_cluster_y.extend(spatial_data_cluster_y0)
					spatial_data_cluster_y1=[int(x)+1 for x in spatial_data_cluster_y0]
					spatial_data_cluster_y2=[int(x)+2 for x in spatial_data_cluster_y0]
					spatial_data_cluster_y3=[int(x)-1 for x in spatial_data_cluster_y0]
					spatial_data_cluster_y4=[int(x)-2 for x in spatial_data_cluster_y0]
					spatial_data_cluster_y.extend(spatial_data_cluster_y1)
					spatial_data_cluster_y.extend(spatial_data_cluster_y2)
					spatial_data_cluster_y.extend(spatial_data_cluster_y3)
					spatial_data_cluster_y.extend(spatial_data_cluster_y4)

				spatial_data_cluster_x = np.array(spatial_data_cluster_x)
				spatial_data_cluster_y = np.array(spatial_data_cluster_y)

				spatial_data_cluster_y_exp=spatial_data_cluster_y
				spatial_data_cluster_x_exp=spatial_data_cluster_x

				# handling boundary conditions:
				spatial_data_cluster_y_exp=np.array([x if x <= image_tmp.shape[1] else image_tmp.shape[1] for x in spatial_data_cluster_y_exp]).astype(int)
				spatial_data_cluster_x_exp=np.array([x if x <= image_tmp.shape[2] else image_tmp.shape[2] for x in spatial_data_cluster_x_exp]).astype(int)

				mask = np.zeros(shape=image_tmp.shape[1:3], dtype='uint8')
				mask[spatial_data_cluster_y_exp-1,spatial_data_cluster_x_exp-1] = 255

				cluster_idx=list(set(metadata['cluster'])).index(list(set(cluster['cluster']))[0])
				if cluster_idx == 0 and self.compress:
					mask_size=(np.array(mask,dtype='uint8').nbytes/1e9)*len(binned_metadata)
					self.compression_factor=np.min([1,np.sqrt(self.target_size/(mask_size+self.raw_size))])
					print('compression factor: '+str(self.compression_factor))

					if self.compression_factor < 0.2:
						print('\n' + '='*60)
						print('⚠️  WARNING: COMPRESSION FACTOR TOO LOW!')
						print('='*60)
						print(f'Compression factor: {self.compression_factor:.4f}')
						print(f'Estimated data loss: {(1-self.compression_factor)*100:.1f}%')
						print(f'Raw size: {self.raw_size:.2f} GB')
						print(f'Mask size: {mask_size:.2f} GB')
						print(f'Total size: {(self.raw_size + mask_size):.2f} GB')
						print(f'Target size: {self.target_size:.2f} GB')
						print('\nThis will result in severe data loss and unusable output!')
						print('AUTOMATICALLY DISABLING COMPRESSION to preserve data.')
						print('\nRecommended solutions for future runs:')
						print(f'  1. Increase target_size to at least {(self.raw_size + mask_size)/5:.1f} GB')
						print('  2. Disable compression (set compress=False)')
						print('  3. Reduce number of mask channels')
						print('='*60 + '\n')
						self.compress = False
						self.compression_factor = 1.0

				cluster_channel=str(cluster_id)+'_mask-expanded'

				if self.compress:
					if self.compression_factor < 0.1:
						print(f'⚠️  Compression factor {self.compression_factor} is dangerously low. Skipping compression for this mask.')
						channels[cluster_channel] = mask
					elif self.filter_img:
						masked_image_preFiltered=ndimage.interpolation.spline_filter1d(mask)
						masked_image_filtered=zoom(masked_image_preFiltered,self.compression_factor)
						channels[cluster_channel]=masked_image_filtered
					else:
						masked_image_nf=zoom(mask,self.compression_factor)
						cluster_channel=str(cluster_id)+'_mask-expanded-nf'
						print('Appending expanded compressed, non-filtered image')
						channels[cluster_channel]=masked_image_nf

				else:
					masked_image = mask
					print('Appending expanded image')
					channels[cluster_channel]=mask
					print(np.sum(image))
			else:
				skipped_clusters += 1

		print(f"\n      Mask generation complete:")
		print(f"      - Valid clusters processed: {valid_clusters}")
		print(f"      - Clusters skipped (≤1 cell): {skipped_clusters}")
		print(f"      - Total mask channels created: {len(channels)}")

		if len(channels) == 0:
			print("      ⚠️  WARNING: No mask channels were created!")
			print("      Possible reasons:")
			print("      - All clusters have ≤1 cell")
			print("      - Spatial annotations are empty or malformed")

		return channels

	@staticmethod
	def _build_ome_xml(channel_names, shape, dtype):
		"""
		Manually construct OME-XML with explicit Channel Name and
		SamplesPerPixel=1 per channel. Used for the >=4GB path so
		bfconvert reads and preserves channel names correctly.
		"""
		import html

		n_channels, size_y, size_x = shape
		dtype_map = {'uint8': 'uint8', 'uint16': 'uint16', 'float32': 'float'}
		ome_dtype = dtype_map.get(str(dtype), 'uint8')

		channel_lines = '\n      '.join(
			'<Channel ID="Channel:0:{i}" Name="{name}" SamplesPerPixel="1"/>'.format(
				i=i, name=html.escape(str(name))
			)
			for i, name in enumerate(channel_names)
		)

		tiffdata_lines = '\n      '.join(
			'<TiffData FirstC="{i}" FirstZ="0" FirstT="0" IFD="{i}" PlaneCount="1"/>'.format(i=i)
			for i in range(n_channels)
		)

		return (
			'<?xml version="1.0" encoding="UTF-8"?>'
			'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"'
			' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
			' xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06'
			' http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">'
			'<Image ID="Image:0" Name="Masquerade">'
			'<Pixels ID="Pixels:0"'
			' Type="{ome_dtype}"'
			' SizeX="{size_x}" SizeY="{size_y}" SizeC="{n_channels}"'
			' SizeZ="1" SizeT="1"'
			' DimensionOrder="XYCZT"'
			' Interleaved="false"'
			' BigEndian="false">'
			'{channel_lines}'
			'{tiffdata_lines}'
			'</Pixels>'
			'</Image>'
			'</OME>'
		).format(
			ome_dtype=ome_dtype,
			size_x=size_x,
			size_y=size_y,
			n_channels=n_channels,
			channel_lines=channel_lines,
			tiffdata_lines=tiffdata_lines,
		)

	@staticmethod
	def write_ome_bigTiff(channels, out, channels_to_keep, tile_size=512, compression=None, pyramid_levels=None, write_mode='ome'):
		"""
		Write a multi-channel ImageJ TIFF readable by QuPath.

		imagej=True and bigtiff=True are mutually exclusive in tifffile —
		bigtiff silently drops the ImageJ metadata block, leaving a plain
		multi-page TIFF that QuPath reads as RGB. Instead we use deflate
		compression to keep the file under the 4GB ImageJ limit. Mask
		channels are ~15% non-zero binary data and compress heavily.
		"""
		import os

		out_dir = os.path.dirname(out)
		if out_dir and not os.path.exists(out_dir):
			print(f'Creating output directory: {out_dir}')
			os.makedirs(out_dir, exist_ok=True)

		new_tiff = np.stack([channels[x] for x in channels_to_keep])
		new_tiff = np.array(new_tiff, dtype='uint8')
		channel_names = channels_to_keep

		uncompressed_gb = new_tiff.nbytes / 1e9
		use_bigtiff = uncompressed_gb >= 4.0

		print(f"Writing ImageJ TIFF for QuPath")
		print(f"  Output: {out}")
		print(f"  Shape: {new_tiff.shape} (C, Y, X)")
		print(f"  Channels: {len(channel_names)}")
		print(f"  dtype: {new_tiff.dtype}")
		print(f"  Uncompressed size: {uncompressed_gb:.2f} GB")
		print(f"  BigTIFF: {use_bigtiff}")
		print(f"  Non-zero: {np.count_nonzero(new_tiff):,} / {new_tiff.size:,} ({100*np.count_nonzero(new_tiff)/new_tiff.size:.3f}%)")

		if np.count_nonzero(new_tiff) == 0:
			print("  ⚠️  CRITICAL WARNING: Output data is completely empty!")

		write_kwargs = {
			'imagej':      True,
			'photometric': 'minisblack',
			'tile':        (512, 512),
			'metadata':    {'Labels': channel_names},
		}
		if use_bigtiff:
			# >= 4GB: BigTIFF with manual OME-XML — bfconvert will pyramidalize
			ome_xml = Masquerade._build_ome_xml(channel_names, new_tiff.shape, new_tiff.dtype)
			page_kwargs = {'photometric': 'minisblack', 'metadata': None}
			if tile_size is not None:
				page_kwargs['tile'] = (tile_size, tile_size)
			with tifffile.TiffWriter(out, bigtiff=True) as tif:
				for i, plane in enumerate(new_tiff):
					tif.write(
						plane,
						description=ome_xml if i == 0 else None,
						**page_kwargs,
					)
		else:
			# < 4GB: same manual OME-XML, no bigtiff.
			# Output must be .ome.tiff so QuPath uses Bio-Formats reader
			# which reads the channel names from the OME-XML.
			ome_xml = Masquerade._build_ome_xml(channel_names, new_tiff.shape, new_tiff.dtype)
			page_kwargs = {'photometric': 'minisblack', 'metadata': None}
			if tile_size is not None:
				page_kwargs['tile'] = (tile_size, tile_size)
			with tifffile.TiffWriter(out) as tif:
				for i, plane in enumerate(new_tiff):
					tif.write(
						plane,
						description=ome_xml if i == 0 else None,
						**page_kwargs,
					)

		written_gb = os.path.getsize(out) / 1e9
		print(f"  ✓ Wrote {len(channel_names)}-channel {'OME-TIFF (will be pyramidalized by bfconvert)' if use_bigtiff else 'ImageJ TIFF'}")
		print(f"  ✓ Written size: {written_gb:.2f} GB")

		with tifffile.TiffFile(out) as tif:
			print(f"  ✓ Pages: {len(tif.pages)}")
			if use_bigtiff:
				print(f"  ✓ Is OME: {tif.is_ome}")
				if tif.ome_metadata:
					from xml.etree import ElementTree as ET
					root = ET.fromstring(tif.ome_metadata)
					ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
					found = root.findall('.//ome:Channel', ns)
					print(f"  ✓ OME channel elements: {len(found)}")
					for ch in found[:3]:
						print(f"    - {ch.get('Name', '(unnamed)')}")
					if len(found) > 3:
						print(f"    ... and {len(found)-3} more")
			else:
				print(f"  ✓ Is ImageJ: {tif.is_imagej}")
				if tif.imagej_metadata:
					labels = tif.imagej_metadata.get('Labels', [])
					print(f"  ✓ Labels in metadata: {len(labels)}")
					if labels:
						print(f"    First: {labels[0]}")
						print(f"    Last:  {labels[-1]}")

		return None

	def get_continuous_channels(self, channels, set_subset_x=None, set_subset_y=None):
		"""
		Updated to handle both OME-TIFF and QPTIFF formats
		"""
		if self.tiff_type is None:
			self.tiff_type = self.detect_tiff_type()

		if self.relevant_markers is not None:
			relevant_markers_ = pd.read_csv(self.relevant_markers)
			markers_use = relevant_markers_['x']
			markers_mod = ["".join(re.split('_',x)) for x in markers_use]
			markers_mod_tmp = [re.sub('_','-',x) for x in markers_use]
			markers_mod.extend(markers_mod_tmp)

		print('pulling input qptiff channels')
		channel_names = []

		with TiffFile(self.image_source) as tif:
			if self.tiff_type == 'ome':
				all_channel_names = self.get_channel_names_ome(tif)
			elif self.tiff_type == 'qptiff':
				all_channel_names = self.get_channel_names_qptiff(tif)
			else:
				all_channel_names = [f"Channel_{i}" for i in range(len(tif.series[0].pages))]

			for page_idx, page in enumerate(tif.series[0].pages):
				page_name = all_channel_names[page_idx] if page_idx < len(all_channel_names) else f"Channel_{page_idx}"

				print(f"{page_idx}\t{page_name}")

				# Handle duplicate channel names
				if page_name in channels:
					page_name = page_name + '_'+str(channel_names.count(page_name)+1)
				channel_names.append(page_name)

				should_include = True
				if self.relevant_markers is not None:
					channel_mod = "".join(re.split("-",page_name))
					should_include = (page_name in markers_use or page_name in markers_mod)

				if should_include:
					page_tmp = page.asarray()

					if self.adjust_coords:
						page_tmp = page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
						page_tmp = page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]

					channels[page_name] = page_tmp

		return channels, channel_names
