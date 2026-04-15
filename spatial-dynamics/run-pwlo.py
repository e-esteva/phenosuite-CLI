from pwlo_es_pt import *
import sys
import pandas as pd

spatial_obj_path = sys.argv[1]
out_dir = sys.argv[2]
label = sys.argv[3]
resolution = float(sys.argv[4])
p1 = float(sys.argv[5])
p2 = float(sys.argv[6])

spatial_obj = pd.read_csv(spatial_obj_path)
pairwise_logOdds(spatial_obj,out_dir,label,draw=False,resolution=resolution,p1=p1,p2=p2,compute_effect_size=False)
