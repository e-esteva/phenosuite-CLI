def pairwise_logOdds(spatial_obj,out_dir,label,draw=False,resolution=0.3774,p1=3,p2=30,compute_effect_size=False):
  import numpy as np
  import pandas as pd
  import scipy
  from scipy import spatial
  from scipy.stats import ks_2samp
  import os
  from itertools import chain
  
  if compute_effect_size:
    c_l=1e3
    downsampling_factor=1/(c_l)
  else:
    c_l=1e4
  
  print(spatial_obj)
  print(compute_effect_size)
  set_sizes=spatial_obj['cluster'].value_counts()
  celltypes=set_sizes.keys()
  
  # resolution=0.377
  # minimum distance:
  #p1=3
  # maximum distance
  #p2=30
  p2=p2+p1
  
  p1_scaled=p1/resolution
  p2_scaled=p2/resolution
  
  compute_ks=compute_effect_size
  #combination_method=['stouffer','fisher'][0]
  
  global_logodds=[]
  global_effect_sizes=[]
  global_probabilities=[]
  for i in range(len(celltypes)):
  
    neighbors_i=[]
  
    ref_=celltypes[i]
    ref_cl=set_sizes[ref_] > c_l
    ref_cl_ratio=set_sizes[ref_]/c_l
    ref_data=spatial_obj[['x','y']][(spatial_obj['cluster']==ref_)]
    
    ref_int_depth=int(np.floor(ref_cl_ratio))
    ref_mod_depth=int(set_sizes[ref_] % c_l)
    
    # wlll contain the Kolmogorov-Smirnov effect_sizes
    blocks_ks=[]
    # number of cells in each block that have interactions with query over all blocks:
    counts_i=[]
    
    for j in range(len(celltypes)):
      query_=celltypes[j]
      print(str(ref_) +'| N= '+str(set_sizes[ref_])+' : '+str(query_)+'| N= '+str(set_sizes[query_]))
      
      query_cl=set_sizes[query_] > c_l
      query_cl_ratio=set_sizes[query_]/c_l
      query_data=spatial_obj[['x','y']][(spatial_obj['cluster']==query_)]
      
      query_int_depth=int(np.floor(query_cl_ratio))
      query_mod_depth=int(set_sizes[query_] % c_l)
      
      
      
      if set_sizes[query_]*set_sizes[ref_] > c_l**2:
       
        if query_int_depth >= 1 and ref_int_depth >= 1:
           
          query_int_blocks=[query_data[int(x*c_l):int(c_l*(x+1))] for x in range(query_int_depth)]
          query_mod_block=query_data[int(query_int_depth*c_l):query_data.shape[0]]
            
          ref_int_blocks=[ref_data[int(x*c_l):int(c_l*(x+1))] for x in range(ref_int_depth)]
          ref_mod_block=ref_data[int(ref_int_depth*c_l):ref_data.shape[0]]
            
            
          # for every ref int block, compare to all query int blocks:
          int_block_neighbors=[]
          mod_block_neighbors_ref_int=[]
          ref_mod_query_int_blocks=[]
          
          ks_blocks=[]
          ks_weights=[]
          
          counts_int_block_neighbors=[]
          counts_mod_block_neighbors_ref_int=[]
          counts_ref_mod_query_int_blocks=[]
          
          for a in range(len(ref_int_blocks)):
            ref_query_blocks=[]
            
            for b in range(len(query_int_blocks)):
              # ref int to query int:
              ref_q_int_dist=spatial.distance.cdist(ref_int_blocks[a],query_int_blocks[b])
              ref_query_blocks.append(ref_q_int_dist)
              
              if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(c_l))
                random_block_dist=spatial.distance.cdist(ref_int_blocks[a],random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor)))
                  
                ref_q_int_dist=np.concatenate(ref_q_int_dist)
                ref_q_int_dist_s=np.random.choice(a=ref_q_int_dist,size=int(np.round(len(ref_q_int_dist)*downsampling_factor)))
                  
                  
                ks_blocks.append(ks_2samp(ref_q_int_dist_s,rb_s,alternative='greater').statistic)
                  
            # now compare ref int block to query mod block:
            ref_int_query_mod_dist=spatial.distance.cdist(ref_int_blocks[a],query_mod_block)
            ref_query_blocks.append(ref_int_query_mod_dist)
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_mod_block.shape[0]))
                random_block_dist=spatial.distance.cdist(ref_int_blocks[a],random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*10)))
                  
                ref_int_query_mod_dist=np.concatenate(ref_int_query_mod_dist)
                ref_int_query_mod_dist_s=np.random.choice(a=ref_int_query_mod_dist,size=int(np.round(len(ref_int_query_mod_dist)*downsampling_factor*10)))
                  
                ks_blocks.append(ks_2samp(ref_int_query_mod_dist_s,rb_s,alternative='greater').statistic)
              
              
            int_block_neighbors_ = [np.sum([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
            int_block_neighbors.append(int_block_neighbors_)
            print(len(int_block_neighbors))
            
            # length: number of blocks. Each entry contains number of ref cells that have interactions
            # produces binary vector across all ref cells indicating presence of a neighbor:
            counts_int_block_neighbors_ = [np.sign([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
            # takes index for ref cells with neighbor:
            counts_int_block_neighbors_clean =[[x for x in range(ref_query_blocks[i].shape[0]) if counts_int_block_neighbors_[i][x] == 1] for i in range(len(ref_query_blocks))]
            # converts lists to tuples for set calculation:
            set_ = [tuple(x) for x in counts_int_block_neighbors_clean]
            # appends the number of uniq ref idx:
            counts_int_block_neighbors.append(len(set(list(chain.from_iterable(set_)))))
            
              
          # now compare ref mod to query int:
          for b in range(len(query_int_blocks)):
            # ref mod to query int:
            ref_mod_q_int_dist=spatial.distance.cdist(ref_mod_block,query_int_blocks[b])
            ref_mod_query_int_blocks.append(ref_mod_q_int_dist)
            
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_int_blocks[b].shape[0]))
                random_block_dist=spatial.distance.cdist(ref_mod_block,random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*10)))
                  
                ref_mod_q_int_dist=np.concatenate(ref_mod_q_int_dist)
                ref_mod_q_int_dist_s=np.random.choice(a=ref_mod_q_int_dist,size=int(np.round(len(ref_mod_q_int_dist)*downsampling_factor*10)))
                  
                ks_blocks.append(ks_2samp(ref_mod_q_int_dist_s,rb_s,alternative='greater').statistic)
                  
          # now compare ref mod block to query mod block:
          ref_mod_query_mod_dist=spatial.distance.cdist(ref_mod_block,query_mod_block)
          ref_mod_query_int_blocks.append(ref_mod_query_mod_dist)
          
          if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_mod_block.shape[0]))
                random_block_dist=spatial.distance.cdist(ref_mod_block,random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb))))
                  
                ref_mod_query_mod_dist=np.concatenate(ref_mod_query_mod_dist)
                ref_mod_query_mod_dist_s=np.random.choice(a=ref_mod_query_mod_dist,size=int(np.round(len(ref_mod_query_mod_dist))))
                  
                ks_blocks.append(ks_2samp(ref_mod_q_int_dist_s,rb_s,alternative='greater').statistic)
              
          mod_block_neighbors_ = [np.sum([np.sum((ref_mod_query_int_blocks[x][z] > p1_scaled) & (ref_mod_query_int_blocks[x][z] < p2_scaled)) for z in range(ref_mod_query_int_blocks[x].shape[0])]) for x in range(len(ref_mod_query_int_blocks))]
          
          #mod_block_neighbors_counts_ = [np.sum(np.sign([np.sum((ref_mod_query_int_blocks[x][z] > p1_scaled) & (ref_mod_query_int_blocks[x][z] < p2_scaled)) for z in range(ref_mod_query_int_blocks[x].shape[0])])) for x in range(len(ref_mod_query_int_blocks))]
          counts_mod_block_neighbors_ = [np.sign([np.sum((ref_mod_query_int_blocks[x][z] > p1_scaled) & (ref_mod_query_int_blocks[x][z] < p2_scaled)) for z in range(ref_mod_query_int_blocks[x].shape[0])]) for x in range(len(ref_mod_query_int_blocks))]
      
          counts_mod_block_neighbors_clean =[[x for x in range(ref_mod_query_int_blocks[i].shape[0]) if counts_mod_block_neighbors_[i][x] == 1] for i in range(len(ref_mod_query_int_blocks))]
          set_ = [tuple(x) for x in counts_mod_block_neighbors_clean]
          
          # int_block_neighbors up till now will have len #ref int block
          # with each entry being of length #query int block +1 (query mod)
          # now: will have len #ref int block +1 (this last entry will have ref mod v all query int + query mod)
          int_block_neighbors.append(mod_block_neighbors_)
          print(len(int_block_neighbors))
          
          counts_int_block_neighbors.append(len(set(list(chain.from_iterable(set_)))))
          
          
          neighbors_i.append(np.sum([np.sum(x) for x in int_block_neighbors]))
          
          counts_i.append(np.sum(counts_int_block_neighbors))
          
          if compute_ks:
            effect_size_combined=np.mean(np.array(ks_blocks))
            
            blocks_ks.append(effect_size_combined)
            print( 'combined effect size: '+str(effect_size_combined))
            
        elif ref_int_depth >= 1:
          ks_blocks=[]
          ks_weights=[]
          
          ref_int_blocks=[ref_data[int(x*c_l):int(c_l*(x+1))] for x in range(ref_int_depth)]
          ref_mod_block=ref_data[int(ref_int_depth*c_l):ref_data.shape[0]]
          
          # for every ref block, compare to query data:
          int_block_neighbors=[]
          mod_block_neighbors_ref_int=[]
          
          ref_mod_query_int_blocks=[]
          ref_query_blocks=[]
          
          counts_int_block_neighbors=[]
          counts_mod_block_neighbors_ref_int=[]
          
          for a in range(len(ref_int_blocks)):
            
            ref_int_query_data_dist=spatial.distance.cdist(ref_int_blocks[a],query_data)
            ref_query_blocks.append(ref_int_query_data_dist)
            
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_data.shape[0]))
                random_block_dist=spatial.distance.cdist(ref_int_blocks[a],random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*10)))
                  
                ref_int_query_data_dist=np.concatenate(ref_int_query_data_dist)
                ref_int_query_data_dist_s=np.random.choice(a=ref_int_query_data_dist,size=int(np.round(len(ref_int_query_data_dist)*downsampling_factor*10)))
                  
                ks_blocks.append(ks_2samp(ref_int_query_data_dist_s,rb_s,alternative='greater').statistic)
            
            
            
          if ref_mod_block.shape[0] > 0:
            # now compare ref mod block to query data:
            ref_mod_query_dist=spatial.distance.cdist(ref_mod_block,query_data)
            ref_query_blocks.append(ref_mod_query_dist)
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_data.shape[0]))
                random_block_dist=spatial.distance.cdist(ref_mod_block,random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*100)))
                  
                ref_mod_query_dist=np.concatenate(ref_mod_query_dist)
                ref_mod_query_dist_s=np.random.choice(a=ref_mod_query_dist,size=int(np.round(len(ref_mod_query_dist)*downsampling_factor*100)))
                  
                ks_blocks.append(ks_2samp(ref_int_query_data_dist_s,rb_s,alternative='greater').statistic)
            
            
          int_block_neighbors_ = [np.sum([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
          int_block_neighbors.append(int_block_neighbors_)    
          
          counts_int_block_neighbors_ = [np.sign([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
          counts_int_block_neighbors_clean =[[x for x in range(ref_query_blocks[i].shape[0]) if counts_int_block_neighbors_[i][x] == 1] for i in range(len(ref_query_blocks))]
          set_ = [tuple(x) for x in counts_int_block_neighbors_clean]
          
          counts_int_block_neighbors.append(len(set(list(chain.from_iterable(set_)))))   
          
          # summing over blocks
          neighbors_i.append(np.sum([np.sum(x) for x in int_block_neighbors]))
          counts_i.append(np.sum(counts_int_block_neighbors))
          
          if compute_ks:
            effect_size_combined=np.mean(np.array(ks_blocks))
            
            blocks_ks.append(effect_size_combined)
            print( 'combined effect size: '+str(effect_size_combined))
            
            
        elif query_int_depth >= 1:
          
          ks_blocks=[]
          ks_weights=[]
          
          query_int_blocks=[query_data[int(x*c_l):int(c_l*(x+1))] for x in range(query_int_depth)]
          query_mod_block=query_data[int(query_int_depth*c_l):query_data.shape[0]]
          
          # for every query block, compare to refdata:
          int_block_neighbors=[]
          mod_block_neighbors_ref_int=[]
          
          ref_mod_query_int_blocks=[]
          ref_query_blocks=[]
          
          counts_int_block_neighbors=[]
          counts_mod_block_neighbors_ref_int=[]
          
          for a in range(len(query_int_blocks)):
            
            query_int_refdata_dist=spatial.distance.cdist(ref_data,query_int_blocks[a])
            ref_query_blocks.append(query_int_refdata_dist)
            
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_int_blocks[a].shape[0]))
                random_block_dist=spatial.distance.cdist(ref_data,random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*10)))
                  
                query_int_refdata_dist=np.concatenate(query_int_refdata_dist)
                query_int_refdata_dist_s=np.random.choice(a=query_int_refdata_dist,size=int(np.round(len(query_int_refdata_dist)*downsampling_factor*10)))
                print(str(ks_2samp(query_int_refdata_dist_s,rb_s,alternative='greater')))
                ks_blocks.append(ks_2samp(query_int_refdata_dist_s,rb_s,alternative='greater').statistic)
          
          if query_mod_block.shape[0] > 0:
            # now compare query mod block to refdata:
            ref_data_query_mod_dist=spatial.distance.cdist(ref_data,query_mod_block)
            ref_query_blocks.append(ref_data_query_mod_dist)
            
            if compute_ks:
                random_block=spatial_obj[['x','y']].sample(int(query_mod_block.shape[0]))
                random_block_dist=spatial.distance.cdist(ref_data,random_block)
                
                
                rb=np.concatenate(random_block_dist)
                rb_s=np.random.choice(a=rb,size=int(np.round(len(rb)*downsampling_factor*100)))
                  
                ref_data_query_mod_dist=np.concatenate(ref_data_query_mod_dist)
                ref_data_query_mod_dist_s=np.random.choice(a=ref_data_query_mod_dist,size=int(np.round(len(ref_data_query_mod_dist)*downsampling_factor*100)))
                  
                print(str(ks_2samp(ref_data_query_mod_dist_s,rb_s,alternative='greater')))
                ks_blocks.append(ks_2samp(ref_data_query_mod_dist_s,rb_s,alternative='greater').statistic)
                
            
          int_block_neighbors_ = [np.sum([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
          int_block_neighbors.append(int_block_neighbors_)    
          
          counts_int_block_neighbors_ = [np.sign([np.sum((ref_query_blocks[x][z] > p1_scaled) & (ref_query_blocks[x][z] < p2_scaled)) for z in range(ref_query_blocks[x].shape[0])]) for x in range(len(ref_query_blocks))]
          counts_int_block_neighbors_clean =[[x for x in range(ref_query_blocks[i].shape[0]) if counts_int_block_neighbors_[i][x] == 1] for i in range(len(ref_query_blocks))]
          set_ = [tuple(x) for x in counts_int_block_neighbors_clean]
          
          counts_int_block_neighbors.append(len(set(list(chain.from_iterable(set_)))))
          
          neighbors_i.append(np.sum([np.sum(x) for x in int_block_neighbors]))
          counts_i.append(np.sum(counts_int_block_neighbors))
          
          if compute_ks:
            effect_size_combined=np.mean(np.array(ks_blocks))
            
            blocks_ks.append(effect_size_combined)
            print( 'combined effect size: '+str(effect_size_combined))
            
            
      else:
          # ref x query dist: each index a row in ref over all query cells:
          c_dist_ij=spatial.distance.cdist(ref_data,query_data)
          neighbors_i.append(np.sum([np.sum((c_dist_ij[x] > p1_scaled) & (c_dist_ij[x] < p2_scaled)) for x in range(c_dist_ij.shape[0])]))
          
          counts_i.append(np.sum(np.sign([np.sum((c_dist_ij[x] > p1_scaled) & (c_dist_ij[x] < p2_scaled)) for x in range(c_dist_ij.shape[0])])))
          
          if compute_ks:
                  random_block=spatial_obj[['x','y']].sample(int(query_data.shape[0]))
                  random_block_dist=spatial.distance.cdist(ref_data,random_block)
                  ks_block=ks_2samp(np.concatenate(c_dist_ij),np.concatenate(random_block_dist),alternative='greater')
                  print('effect-size: '+str(ks_block.statistic))
                  blocks_ks.append(ks_block.statistic)       
        
      
      
    log_odds_i=np.array([np.log( (neighbors_i[x]/np.sum(np.array(neighbors_i))) / (set_sizes[x]/spatial_obj.shape[0])) for x in range(len(neighbors_i))]) 
    print('log_odds_i: '+str(log_odds_i))
    global_logodds.append(log_odds_i)
    
    probabilities_i = np.array(counts_i)/set_sizes[ref_]
    
    global_probabilities.append(probabilities_i)
    
    if compute_ks:
      global_effect_sizes.append(blocks_ks)
      
  ###################################################  
  global_logodds_df=pd.DataFrame(global_logodds).T
  global_logodds_df.columns=celltypes
  global_logodds_df.index=celltypes
      
  global_logodds_df.to_csv(str(out_dir)+'/'+str(label)+'-logOdds_matrix.csv')
  
  # probabilities:
  global_probabilities_df=pd.DataFrame(global_probabilities).T
  global_probabilities_df.columns=celltypes
  global_probabilities_df.index=celltypes
      
  global_probabilities_df.to_csv(str(out_dir)+'/'+str(label)+'-probabilities_matrix.csv')
  
  ## Effect-sizes:
  if compute_ks:
      global_effect_sizes_df=pd.DataFrame(global_effect_sizes).T
      global_effect_sizes_df.columns=celltypes
      global_effect_sizes_df.index=celltypes
      
      global_effect_sizes_df.to_csv(str(out_dir)+'/'+str(label)+'-KS-effect_sizes_matrix.csv')
   
  return global_logodds_df
