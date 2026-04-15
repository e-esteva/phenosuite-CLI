import sys
from getNeighborhoods import *
from optimize_neighborhoods import *


spatial_obj = sys.argv[1]
out_dir = sys.argv[2]
label = sys.argv[3]
target_celltypes = sys.argv[4]

neighborhood_metadata =  getNeighborhoods(spatial_obj=spatial_obj, out_dir=out_dir, label=label, target_celltypes=target_celltypes, target_clustering='cluster')

optimize_neighborhood(neighborhood_metadata['neighborhoods'], out_dir, thresholds=np.arange(0.5, 0.95, 0.05), make_gif=True, min_cells=100)
