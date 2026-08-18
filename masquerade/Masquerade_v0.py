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
	
	def pullExt(self, file):
		return re.split('[.]',file)[-1]
	
	def PreProcessImage(self):
		# import spatial metadata:
		spatial_metadata=pd.read_csv(self.spatial_anno)
		
		# import qptiff:
		image = skimage.io.imread(fname=str(self.image_source),plugin='tifffile')
		# load in cluster spatial coords & subset image to range captured by coords:
		
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
		channel_names=[]
		with TiffFile(self.image_source) as tif:
			for page in tif.series[0].pages:
				page_name_element = ElementTree.fromstring(page.description).find('Biomarker')
				if page_name_element is not None:
					page_name = page_name_element.text
					page_name = page_name.replace(' ', '-')
					print(page.index, '\t', page_name)
					
					if page_name in channels:
						page_name = page_name + '_'+str(channel_names.count(page_name)+1)
						channel_names.append(page_name)
					
					channel_mod="".join(re.split("-",page_name))
					if self.relevant_markers is not None:	
						relevant_markers_ = pd.read_csv(self.relevant_markers)
						markers_use=relevant_markers_['x']
						markers_mod=["".join(re.split('_',x)) for x in markers_use]
						markers_mod_tmp=[re.sub('_','-',x) for x in markers_use]
						markers_mod.extend(markers_mod_tmp)
						
						if page_name in markers_use or page_name in markers_mod:
							page_tmp=page.asarray()
							if self.adjust_coords:
								subset_x=np.array(spatial_metadata['x']).astype(int)
								set_subset_x=np.array([x for x in set(subset_x)]).astype(int)
								
								subset_y=np.array(spatial_metadata['y']).astype(int)
								set_subset_y=np.array([x for x in set(subset_y)]).astype(int)
								
								page_tmp=page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
								page_tmp=page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]
							
							page_tmp=zoom(page_tmp,self.compression_factor)
							channels[page_name] = page_tmp
					else:
						page_tmp=page.asarray()
						if self.adjust_coords:
							subset_x=np.array(spatial_metadata['x']).astype(int)
							set_subset_x=np.array([x for x in set(subset_x)]).astype(int)
							
							subset_y=np.array(spatial_metadata['y']).astype(int)
							set_subset_y=np.array([x for x in set(subset_y)]).astype(int)
							
							page_tmp=page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
							page_tmp=page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]
							page_tmp=zoom(page_tmp,self.compression_factor)
							channels[page_name] = page_tmp
						
						else:
							page_tmp=zoom(page_tmp,self.compression_factor)
							channels[page_name] = page_tmp
		
		return channels,channel_names
	
	@staticmethod
	def writeMaskTiff(channels, outPath):
		num_channels=len(channels.keys())
		print('Markers retained: ')
		print(channels.keys())
		
		# stacking all channels:
		new_tiff=np.stack([channels[x] for x in channels.keys()])
		new_tiff=np.array(new_tiff,dtype='uint8')
		
		print("writing new tiff")
		channels=list(channels.keys())
		ijmetadata={"Labels":channels}
		
		# write new tiff:
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
		# Generate angles
		theta = np.linspace(0, 2*np.pi, num_points)
		
		# Generate points on circle
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
		# dynamic design module:
		binned_metadata=[metadata[metadata['cluster']==x] for x in list(set(metadata['cluster']))]
		
		masks = {}
		masks_preFiltered = {}
		channels = {}
		
		# circle filled: presentation; empty: phenotyping 
		for cluster in tqdm.tqdm(binned_metadata):
			if cluster.shape[0] > 1:
				cluster_id=cluster['cluster'][cluster['cluster'].keys()[0]]
				
				if self.adjust_coords:
					# adjusting coordinates for subset space:
					cluster_x=cluster['x']-np.min(set_subset_x)
					cluster_y=cluster['y']-np.min(set_subset_y)
				else:
					cluster_x=cluster['x']
					cluster_y=cluster['y']
				
				print('Processing cluster: '+str(cluster_id)+'\n')
				spatial_data_cluster_x=np.array([int(x) for x in cluster_x])
				
				spatial_data_cluster_y=np.array([int(x) for x in cluster_y])        
				
				image_tmp=image
				
				cluster_channel=str(cluster_id)+'_mask-original'
				
				print('Expanding image')
				
				spatial_data_cluster_x0=[int(x) for x in cluster_x]
				spatial_data_cluster_y0=[int(x) for x in cluster_y]
				
				if self.filled:
					circle_coords = [self.generate_filled_circle(cx=spatial_data_cluster_x0[x], cy=spatial_data_cluster_y0[x], r=self.radius) for x in range(len(spatial_data_cluster_x0))]
					
					# x coords: first element in tuple (a list of tuples for each centroid):
					spatial_data_cluster = [y for x in circle_coords for y in x]
					
					spatial_data_cluster_x = [list(x)[0] for x in spatial_data_cluster]
					
					# y coords: 2nd element in tuple:
					spatial_data_cluster_y = [list(x)[1] for x in spatial_data_cluster]
				
				else:
					circle_coords = [self.generate_circle(x_center=spatial_data_cluster_x0[x], y_center=spatial_data_cluster_y0[x], radius=self.radius, num_points=self.num_points) for x in range(len(spatial_data_cluster_x0))]
					# x coords: 1st array in tuple:
					spatial_data_cluster_x_ = [x[0] for x in circle_coords]
					
					spatial_data_cluster_x = [y for x in spatial_data_cluster_x_ for y in x]
					
					# y coords: 2nd array in tuple:
					spatial_data_cluster_y_ = [x[1] for x in circle_coords]
					
					spatial_data_cluster_y = [y for x in spatial_data_cluster_y_ for y in x]
					
					# now add back centroid + 2 in all directions (cross)
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
				
				# handling boundary conditions (literally):
				spatial_data_cluster_y_exp=np.array([x if x <= image_tmp.shape[1] else image_tmp.shape[1] for x in spatial_data_cluster_y_exp]).astype(int)
				
				spatial_data_cluster_x_exp=np.array([x if x <= image_tmp.shape[2] else image_tmp.shape[2] for x in spatial_data_cluster_x_exp]).astype(int)
				
				# new proposal:
				mask = np.zeros(shape=image_tmp.shape[1:3])
				mask[spatial_data_cluster_y_exp-1,spatial_data_cluster_x_exp-1]=1.0
				
				# all mask channels are the same size so:
				cluster_idx=list(set(metadata['cluster'])).index(list(set(cluster['cluster']))[0])
				if cluster_idx == 0 and self.compress:
					mask_size=(np.array(mask,dtype='uint8').nbytes/1e9)*len(binned_metadata)
					self.compression_factor=np.min([1,np.sqrt(self.target_size/(mask_size+self.raw_size))])
					print('compression factor: '+str(self.compression_factor))
				
				cluster_channel=str(cluster_id)+'_mask-expanded'
				
				if self.compress:
					if self.filter_img:
						# compressed filtered
						masked_image_preFiltered=ndimage.interpolation.spline_filter1d(mask)
						masked_image_filtered=zoom(masked_image_preFiltered,self.compression_factor)
						channels[cluster_channel]=masked_image_filtered
					else:
						# compressed, no filtering
						masked_image_nf=zoom(mask,self.compression_factor)
						cluster_channel=str(cluster_id)+'_mask-expanded-nf'
						print('Appending expanded compressed, non-filtered image')
						channels[cluster_channel]=masked_image_nf
				
				# no compression, filtered:
				else:
					masked_image = mask
				
					print('Appending expanded image')
				
					channels[cluster_channel]=mask
				
					print(np.sum(image))
		
		return channels
	
	@staticmethod
	def write_ome_bigTiff(channels, out, channels_to_keep):

		# stacking all channels:
		new_tiff=np.stack([channels[x] for x in channels_to_keep])
		new_tiff=np.array(new_tiff,dtype='uint8')
		
		print("writing new tiff")
		channel_names=channels_to_keep
		ijmetadata={"Labels":channel_names}
		
		labs=[x for x in ijmetadata.values()]
		print(labs[0])
		
		tifffile.imwrite(out, new_tiff, metadata={"Channel": {"Name": labs[0]}}, ome=True, bigtiff=True)
		
		return None
	
	def get_continuous_channels(self,channels,set_subset_x=None,set_subset_y=None):
		if self.relevant_markers is not None:
			relevant_markers_ = pd.read_csv(self.relevant_markers)
			markers_use=relevant_markers_['x']
			markers_mod=["".join(re.split('_',x)) for x in markers_use]
			markers_mod_tmp=[re.sub('_','-',x) for x in markers_use]
			markers_mod.extend(markers_mod_tmp)

		print('pulling input qptiff channels')
		# loading in qptiff channels and adding mask as last channel:
		#channels = {}
		channel_names=[]
		with TiffFile(self.image_source) as tif:
			for page in tif.series[0].pages:
				page_name_element = ElementTree.fromstring(page.description).find('Biomarker')
				if page_name_element is not None:
					page_name = page_name_element.text
					page_name = page_name.replace(' ', '-')
					print(page.index, '\t', page_name)
					#channel_names.append(page_name)

					#if page_name in channels: raise Exception('multiple pages/channels with the same name')
					if page_name in channels: page_name = page_name + '_'+str(channel_names.count(page_name)+1)
					channel_names.append(page_name)

					channel_mod="".join(re.split("-",page_name))
					if self.relevant_markers is not None:	
						if page_name in markers_use or page_name in markers_mod:
							page_tmp=page.asarray()
							if self.adjust_coords:
								page_tmp=page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
								page_tmp=page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]
				
					
							channels[page_name] = page_tmp
					else:
						page_tmp=page.asarray()
						if self.adjust_coords:
							page_tmp=page_tmp[np.min(set_subset_y):np.max(set_subset_y),:]
							page_tmp=page_tmp[:,np.min(set_subset_x):np.max(set_subset_x)]

					
							channels[page_name] = page_tmp
				
						else:
					
							channels[page_name] = page_tmp



		return channels,channel_names


