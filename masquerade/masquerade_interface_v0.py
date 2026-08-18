#from masquerade_utils import *
#from masquerade_utils import PreProcessImage
#from masquerade_utils import get_circle_masks
#from masquerade_utils import get_continuous_channels
#from masquerade_utils import write_ome_bigTiff
import sys
import re
import numpy as np

# ${image} ${spatial_metadata} ${outPath} ${relevant_markers} ${adjust_coords} ${compress} ${radius} ${preFilter_masks} ${target_size}
image_source = sys.argv[1]
metadata = sys.argv[2]
out_path = sys.argv[3]
marker_metadata = sys.argv[4]
adjust_coords = sys.argv[5]
compress = sys.argv[6]
radius = sys.argv[7]
filled = sys.argv[8]
num_points = sys.argv[9]
preFilter_masks = sys.argv[10]
target_size = sys.argv[11]


#metadata='/gpfs/data/abl/tric/segmentation/data/PediatricFollicularLymphoma/out-phenomenalist/20230630_PFL_9457_B2_raw_data.qptiff_object_Data_leftpanel/mask-inputs/20230630_PFL_9457_B2_raw_data.qptiff_object_Data_leftpanel-spatial-anno-res7.csv'

#spatial_metadata=pd.read_csv(metadata)

#image_source='/gpfs/data/abl/tric/segmentation/data/PediatricFollicularLymphoma/20230630_PFL_9457_B2_raw_data.qptiff'

#img,raw_size,subset_x,subset_y = PreProcessImage(image_source=image_source, spatial_metadata=spatial_metadata)

#channels = get_circle_masks(image=img,set_subset_x=subset_x,set_subset_y=subset_y,metadata=spatial_metadata)

#marker_metadata='/gpfs/data/abl/tric/segmentation/data/PediatricFollicularLymphoma/out-phenomenalist/20230630_PFL_9457_B2_raw_data.qptiff_object_Data_leftpanel/mask-inputs/20230630_PFL_9457_B2_raw_data.qptiff_object_Data_leftpanel-marker-metadata.csv'

#channels,channel_names = get_continuous_channels(image_source=image_source,channels=channels,set_subset_x=subset_x,set_subset_y=subset_y,relevant_markers=marker_metadata)

#write_ome_bigTiff(channels=channels,out='/gpfs/data/abl/tric/segmentation/data/PediatricFollicularLymphoma/ome_bigtiff_test.tiff',channels_to_keep=list(channels.keys()))

## class implementation:
from Masquerade import Masquerade

masquerade=Masquerade()
if adjust_coords == 'True':
	adjust_coords = True
else:
	adjust_coords = False

if marker_metadata == 'None':
	marker_metadata = None

if compress == 'True':
	compress = True
else:
	compress = False

if preFilter_masks == 'True':
	preFilter_masks = True
else:
	preFilter_masks = False

if filled == 'True':
	filled = True
else:
	filled = False

print('image source --> '+ image_source)
print('spatial annos --> '+ metadata)
print('out path --> '+ out_path)
print('marker metadata --> '+ str(marker_metadata))
print('adjust coords: '+ str(adjust_coords))
print('compress: '+str(compress))
print('radius: '+str(radius))
print('fill circle: '+str(filled))
print('num points: '+str(num_points))
print('pre-filter masks: '+str(preFilter_masks))
print('target size --> '+str(target_size))


masquerade.relevant_markers=marker_metadata
masquerade.spatial_anno=metadata
masquerade.image_source=image_source
masquerade.compress=compress
masquerade.target_size=float(target_size)
masquerade.filter_img=preFilter_masks
masquerade.radius=int(radius)
masquerade.filled=filled
masquerade.num_points=int(num_points)

print('raw size: '+str(masquerade.raw_size))

image,set_subset_x,set_subset_y,spatial_metadata = masquerade.PreProcessImage()

print('raw size: '+str(masquerade.raw_size))

print('compression factor: '+str(masquerade.compression_factor))

channels = masquerade.get_circle_masks(image=image,metadata=spatial_metadata,set_subset_x=set_subset_x,set_subset_y=set_subset_y)


if compress:

	print('compressing masks to '+str(masquerade.target_size)+'GB')
	print('compression factor: '+str(masquerade.compression_factor))

	channels,channel_names = masquerade.compress_marker_channels(channels,spatial_metadata=spatial_metadata)
	print(channel_names)
	for k in channels.keys():
		print(channels[k].shape)

	out_path_p=re.split('[.]',out_path)
	out_path_ = f"{out_path_p[0]}_{np.round(masquerade.compression_factor,3)}_.{out_path_p[1]}"

else:
	channels,channel_names = masquerade.get_continuous_channels(channels,set_subset_x=set_subset_x,set_subset_y=set_subset_y)

masquerade.write_ome_bigTiff(channels=channels,out=out_path,channels_to_keep=list(channels.keys()))
