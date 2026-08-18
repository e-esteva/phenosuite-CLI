import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsRegressor
import re
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy import stats
import logging
from functools import reduce
from itertools import combinations
import math
from n_simplex_neighborhoods import *

def getNeighborhoods(spatial_obj, out_dir, label, target_celltypes, target_clustering, n_wise_logOdds):

    neighborhood_metadata=n_wise_logOdds(spatial_obj = spatial_obj,out_dir = out_dir,label = label,target_celltypes =target_celltypes)

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    
    # Sort celltypes by frequency
    # First, create a frequency table similar to R's table function
    cluster_counts = spatial_obj['cluster'].value_counts()
    target_celltype_counts = {ct: cluster_counts.get(ct, 0) for ct in target_celltypes}
    celltypes_sorted = sorted(target_celltype_counts.keys(), key=lambda x: target_celltype_counts[x], reverse=True)
    
    # Get the target clustering
    res = target_clustering
    
    # Sort clusters by frequency
    clusters = spatial_obj[res].value_counts().sort_values(descending=True)
    
    # Function to process each permutation
    def process_permutation(cluster_name):
        logging.info(f'dropping cluster {cluster_name}')
        # Drop cluster x
        mask = spatial_obj[res] != cluster_name
        spatial_obj_tmp = spatial_obj[mask].copy()
        
        # Check if all target celltypes are present
        cluster_tmp_counts = spatial_obj_tmp['cluster'].value_counts()
        target_celltype_counts_tmp = {ct: cluster_tmp_counts.get(ct, 0) for ct in target_celltypes}
        
        if all(count > 0 for count in target_celltype_counts_tmp.values()):
            # All target celltypes are present
            metadata_tmp = n_wise_logOdds(
                spatial_obj=spatial_obj_tmp,
                out_dir=out_dir,
                label=f'{label}-{cluster_name}',
                target_celltypes=target_celltypes
            )
        else:
            # Not all target celltypes are present, calculate theoretical probability
            total_cells = len(spatial_obj_tmp)
            target_probabilities = {ct: cluster_tmp_counts.get(ct, 0) / total_cells for ct in target_celltypes}
            
            # Calculate combinations
            n = len(target_celltypes)
            combs = list(combinations(range(n), n-1))
            
            # Calculate probability for each combination
            combination_probabilities = 0
            for comb in combs:
                prob_values = [list(target_probabilities.values())[i] for i in comb]
                if not any(np.isnan(prob_values)):
                    combination_probabilities += np.prod(prob_values)
            
            # Calculate theoretical probability
            theoretical_prob = combination_probabilities / sum(filter(lambda x: not np.isnan(x), target_probabilities.values()))
            metadata_tmp = [None, math.log(1/math.exp(theoretical_prob))]
        
        return metadata_tmp
    
    # Process all permutations
    permutations = []
    for cluster_name in clusters.index:
        permutations.append(process_permutation(cluster_name))
    
    # Extract log_exp_ratios
    log_exp_ratios = [perm[1] for perm in permutations]
    
    # Assuming neighborhood_metadata is a result from previous computation not shown in the code
    # For the refactoring, I'll reference it as an input from the calling function
    deltas = neighborhood_metadata[1] - np.array(log_exp_ratios)
    delta_z = (deltas - np.mean(deltas)) / np.std(deltas)
    
    # Assign deltas to clusters
    delta_scaffold_z = np.zeros(len(spatial_obj))
    for i, cluster_name in enumerate(clusters.index):
        mask = spatial_obj[res] == cluster_name
        delta_scaffold_z[mask] = delta_z[i]
    
    # Create coordinates DataFrame
    coords = pd.DataFrame({
        'x': spatial_obj['x'],
        'y': -spatial_obj['y'],  # Note the sign flip
        'delta_z': delta_scaffold_z
    })
    
    # KNN smoothing
    knn = KNeighborsRegressor(n_neighbors=10)
    knn.fit(coords[['x', 'y']], coords['delta_z'])
    coords['delta_z_smooth'] = knn.predict(coords[['x', 'y']])
    
    # Create plot 1
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    cmap = LinearSegmentedColormap.from_list("custom", ["blue", "white", "red"])
    scatter1 = ax1.scatter(coords['x'], coords['y'], c=coords['delta_z'], cmap=cmap)
    plt.colorbar(scatter1, label="Custom Scale")
    ax1.set_aspect(1 / (max(abs(coords['y'])) / max(coords['x'])))
    ax1.set_title('Delta Z Distribution')
    p1 = fig1
    
    # Create targets column
    coords['targets'] = 'else'
    for celltype in target_celltypes:
        # Use regex to handle the special characters like '+'
        pattern = re.compile(re.escape(celltype))
        mask = spatial_obj['cluster'].apply(lambda x: bool(pattern.search(str(x))))
        coords.loc[mask, 'targets'] = celltype
    
    # Create colors for the target celltypes
    # Using a custom color function similar to ggsci::pal_cosmic()
    def pal_cosmic(n):
        # Simple implementation returning n colors
        cosmic_colors = sns.color_palette("husl", n)
        return cosmic_colors
    
    colors = {}
    for i, celltype in enumerate(target_celltypes):
        colors[celltype] = pal_cosmic(len(target_celltypes))[i]
    
    # Print sanity checks
    print(colors)
    print(celltypes_sorted)
    print(coords.head())
    print(coords['targets'].value_counts())
    
    # Create overlay points
    overlay_points_list = []
    for celltype in celltypes_sorted:
        # Use regex to handle the special characters like '+'
        pattern = re.compile(re.escape(celltype))
        mask = coords['targets'].apply(lambda x: bool(pattern.search(str(x))))
        
        overlay_points_tmp = pd.DataFrame({
            'x': coords.loc[mask, 'x'],
            'y': coords.loc[mask, 'y'],
            'color': [colors[celltype]] * sum(mask)
        })
        overlay_points_list.append(overlay_points_tmp)
    
    overlay_points = pd.concat(overlay_points_list)
    
    # Create plot 2
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    scatter2 = ax2.scatter(coords['x'], coords['y'], c=coords['delta_z'], cmap=cmap, s=15, alpha=0.7)
    
    # Add overlay points
    for celltype in celltypes_sorted:
        mask = overlay_points['color'] == colors[celltype]
        if sum(mask) > 0:
            ax2.scatter(
                overlay_points.loc[mask, 'x'], 
                overlay_points.loc[mask, 'y'],
                color=colors[celltype],
                s=30
            )
    
    plt.colorbar(scatter2, label="Custom Scale")
    ax2.set_aspect(1 / (max(abs(coords['y'])) / max(coords['x'])))
    ax2.set_title('Delta Z with Cell Type Overlay')
    p2 = fig2
    
    # Create combined plot (vertical stack)
    fig3, (ax3_top, ax3_bottom) = plt.subplots(2, 1, figsize=(10, 16))
    
    # Recreate p1 on top subplot
    scatter3_top = ax3_top.scatter(coords['x'], coords['y'], c=coords['delta_z'], cmap=cmap)
    plt.colorbar(scatter3_top, ax=ax3_top, label="Custom Scale")
    ax3_top.set_aspect(1 / (max(abs(coords['y'])) / max(coords['x'])))
    ax3_top.set_title('Delta Z Distribution')
    
    # Recreate p2 on bottom subplot
    scatter3_bottom = ax3_bottom.scatter(coords['x'], coords['y'], c=coords['delta_z'], cmap=cmap, s=15, alpha=0.7)
    
    # Add overlay points
    for celltype in celltypes_sorted:
        mask = overlay_points['color'] == colors[celltype]
        if sum(mask) > 0:
            ax3_bottom.scatter(
                overlay_points.loc[mask, 'x'], 
                overlay_points.loc[mask, 'y'],
                color=colors[celltype],
                s=30
            )
    
    plt.colorbar(scatter3_bottom, ax=ax3_bottom, label="Custom Scale")
    ax3_bottom.set_aspect(1 / (max(abs(coords['y'])) / max(coords['x'])))
    ax3_bottom.set_title('Delta Z with Cell Type Overlay')
    
    combined_plot = fig3
    
    # Return results as a dictionary
    return {
        'neighborhood_metadata': neighborhood_metadata,
        'plots': [p1, p2, combined_plot],
        'permutation_metadata': {
            'log_exp_ratios': log_exp_ratios,
            'deltas': deltas
        },
        'neighborhoods': coords
    }
