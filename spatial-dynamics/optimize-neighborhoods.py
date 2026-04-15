import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import logging
import math
import re
import shutil
from PIL import Image
import glob
from sklearn.neighbors import KNeighborsRegressor

def optimize_neighborhood(neighborhoods, out_dir, thresholds=np.arange(0.5, 0.96, 0.05), make_gif=False, min_cells=100):
    """
    Python implementation of the optimize_neighborhood function from R.
    
    Parameters:
    -----------
    neighborhoods : pandas.DataFrame
        DataFrame containing spatial neighborhood data
    out_dir : str
        Output directory for saving results
    thresholds : numpy.ndarray
        Thresholds for delta_z_smooth filtering
    make_gif : bool
        Whether to create animated GIFs of the refinement process
    min_cells : int
        Minimum number of cells required for a valid neighborhood
    
    Returns:
    --------
    pandas.DataFrame
        Optimally refined neighborhood
    """
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    # Create output directories
    os.makedirs(f"{out_dir}/refinement", exist_ok=True)
    os.makedirs(f"{out_dir}/refinement-objective", exist_ok=True)
    
    # Get target cell types (excluding 'else')
    target_celltypes = neighborhoods['targets'].unique()
    target_celltypes = target_celltypes[target_celltypes != 'else']
    
    # Initialize lists for tracking metrics
    delta_z_spreading = []
    target_densities = []
    thresholded_neighborhoods = []
    objective_values = []
    
    # Process each threshold
    for threshold in thresholds:
        # Filter neighborhoods based on delta_z_smooth threshold
        thresh_value = np.quantile(neighborhoods['delta_z_smooth'], threshold)
        tmp_circuit_region = neighborhoods[neighborhoods['delta_z_smooth'] > thresh_value].copy()
        
        # Check if the filtered region meets minimum requirements
        target_counts = tmp_circuit_region['targets'].value_counts()
        target_present = all(ct in target_counts.index for ct in target_celltypes)
        
        if len(tmp_circuit_region) >= min_cells and target_present:
            logging.info(f'N = {len(tmp_circuit_region)}')
            min_target_count = min([target_counts.get(ct, 0) for ct in target_celltypes])
            logging.info(f'Min target N = {min_target_count}')
            
            # Sort cell types by frequency
            neighborhood_counts = neighborhoods['targets'].value_counts()
            celltypes_sorted = sorted(
                target_celltypes, 
                key=lambda x: neighborhood_counts.get(x, 0), 
                reverse=True
            )
            
            # Create color palette similar to ggsci::pal_igv()
            def pal_igv(n):
                return sns.color_palette("Set1", n)
            
            colors = {}
            for i, ct in enumerate(target_celltypes):
                colors[ct] = pal_igv(len(target_celltypes))[i]
            
            # Create overlay points for visualization
            overlay_points_list = []
            for celltype in celltypes_sorted:
                # Use regex to handle special characters like '+'
                pattern = re.compile(re.escape(celltype))
                mask = tmp_circuit_region['targets'].apply(lambda x: bool(pattern.search(str(x))))
                
                if mask.sum() > 0:
                    overlay_points_tmp = pd.DataFrame({
                        'x': tmp_circuit_region.loc[mask, 'x'],
                        'y': tmp_circuit_region.loc[mask, 'y'],
                        'color': [colors[celltype]] * mask.sum()
                    })
                    overlay_points_list.append(overlay_points_tmp)
            
            if overlay_points_list:
                overlay_points = pd.concat(overlay_points_list, ignore_index=True)
                
                # Create plot with cell type overlays
                fig, ax = plt.subplots(figsize=(10, 8))
                ax.scatter(tmp_circuit_region['x'], tmp_circuit_region['y'], s=3, alpha=0.2)
                
                for celltype in celltypes_sorted:
                    pattern = re.compile(re.escape(celltype))
                    mask = overlay_points['color'] == colors[celltype]
                    if mask.sum() > 0:
                        ax.scatter(
                            overlay_points.loc[mask, 'x'],
                            overlay_points.loc[mask, 'y'],
                            color=colors[celltype],
                            s=30,
                            marker='o'
                        )
                
                ax.set_aspect(1 / (max(abs(tmp_circuit_region['y'])) / max(tmp_circuit_region['x'])))
                ax.set_title(f'delta_z_smooth > {round(thresh_value, 3)}')
                plt.tight_layout()
                
                # Save the plot
                plt.savefig(f"{out_dir}/refinement/{threshold}.png")
                plt.close()
                
                # Calculate metrics
                target_celltype_density = sum(tmp_circuit_region['targets'].isin(target_celltypes)) / len(neighborhoods)
                target_densities.append(target_celltype_density)
                
                # Calculate delta_z spreading
                delta_z_q25 = np.quantile(tmp_circuit_region['delta_z'], 0.25)
                delta_z_q50 = np.quantile(tmp_circuit_region['delta_z'], 0.50)
                delta_z_q75 = np.quantile(tmp_circuit_region['delta_z'], 0.75)
                
                delta_z_smooth_q25 = np.quantile(tmp_circuit_region['delta_z_smooth'], 0.25)
                delta_z_smooth_q50 = np.quantile(tmp_circuit_region['delta_z_smooth'], 0.50)
                delta_z_smooth_q75 = np.quantile(tmp_circuit_region['delta_z_smooth'], 0.75)
                
                spread_value = (
                    math.exp(0.5 * delta_z_q50 + 0.25 * (delta_z_q25 + delta_z_q75)) * 
                    math.exp(0.5 * delta_z_smooth_q50 + 0.25 * (delta_z_smooth_q25 + delta_z_smooth_q75))
                )
                delta_z_spreading.append(spread_value)
                
                # Calculate objective function
                objective_tmp = target_densities[-1] * math.log(delta_z_spreading[-1])
                objective_values.append(objective_tmp)
                
                # Plot objective function
                plt.figure(figsize=(6, 4))
                plt.plot(objective_values, 'b-')
                plt.title('Objective Function')
                plt.savefig(f"{out_dir}/refinement-objective/{threshold}.png")
                plt.close()
                
                # Store the neighborhood
                thresholded_neighborhoods.append(tmp_circuit_region)
    
    # Plot metrics
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(delta_z_spreading, 'b-')
    ax1.set_title('Delta Z Spreading')
    
    ax2.plot(target_densities, 'r-')
    ax2.set_title('Target Densities')
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/metrics.png")
    plt.close()
    
    # Get optimal threshold
    objective_global = objective_values
    max_idx = np.argmax(objective_global)
    optimal_threshold = thresholds[max_idx]
    optimally_refined_neighborhood = thresholded_neighborhoods[max_idx]
    
    logging.info(f'Optimal neighborhood found at idx = {max_idx} | threshold: {optimal_threshold}')
    
    # Plot final objective function
    plt.figure(figsize=(8, 6))
    plt.plot(objective_global, 'b-')
    plt.axhline(y=max(objective_global), color='r', linestyle='-')
    plt.axvline(x=max_idx, color='r', linestyle='--')
    plt.title('Final Objective Function')
    plt.savefig(f"{out_dir}/final-objective.png")
    plt.close()
    
    # Function to create color scale similar to R's create_min_zero_max_scale
    def create_min_zero_max_scale(data_values, name="Cluster-wise ES Z"):
        min_val = min(data_values)
        max_val = max(data_values)
        
        # Define colors
        colors = ["blue", "white", "red"]
        
        if min_val >= 0:
            # All positive values - white to red
            cmap = LinearSegmentedColormap.from_list("custom", ["white", "red"])
            norm = plt.Normalize(0, max_val)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            return sm, cmap, norm
        elif max_val <= 0:
            # All negative values - blue to white
            cmap = LinearSegmentedColormap.from_list("custom", ["blue", "white"])
            norm = plt.Normalize(min_val, 0)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            return sm, cmap, norm
        else:
            # Values span across zero
            cmap = LinearSegmentedColormap.from_list("custom", colors)
            norm = plt.Normalize(min_val, max_val)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            return sm, cmap, norm
    
    # Create and save the optimal neighborhood plot
    sm, cmap, norm = create_min_zero_max_scale(optimally_refined_neighborhood['delta_z'])
    
    plt.figure(figsize=(20, 20))
    plt.scatter(
        optimally_refined_neighborhood['x'],
        optimally_refined_neighborhood['y'],
        c=optimally_refined_neighborhood['delta_z'],
        cmap=cmap,
        norm=norm
    )
    plt.colorbar(sm, label="Cluster-wise ES Z")
    plt.gca().set_aspect(1 / (max(abs(optimally_refined_neighborhood['y'])) / max(optimally_refined_neighborhood['x'])))
    plt.title('Optimal Neighborhood')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/optimal_neighborhood.pdf")
    plt.savefig(f"{out_dir}/optimal_neighborhood.png")
    plt.close()
    
    # Create and save the optimal neighborhood plot with cell type overlays
    overlay_points_list = []
    for celltype in celltypes_sorted:
        pattern = re.compile(re.escape(celltype))
        mask = optimally_refined_neighborhood['targets'].apply(lambda x: bool(pattern.search(str(x))))
        
        if mask.sum() > 0:
            overlay_points_tmp = pd.DataFrame({
                'x': optimally_refined_neighborhood.loc[mask, 'x'],
                'y': optimally_refined_neighborhood.loc[mask, 'y'],
                'color': [colors[celltype]] * mask.sum(),
                'celltype': optimally_refined_neighborhood.loc[mask, 'targets']
            })
            overlay_points_list.append(overlay_points_tmp)
    
    if overlay_points_list:
        overlay_points = pd.concat(overlay_points_list, ignore_index=True)
        
        plt.figure(figsize=(20, 20))
        plt.scatter(
            optimally_refined_neighborhood['x'],
            optimally_refined_neighborhood['y'],
            s=3,
            alpha=0.2
        )
        
        for celltype in celltypes_sorted:
            mask = overlay_points['celltype'] == celltype
            if mask.sum() > 0:
                plt.scatter(
                    overlay_points.loc[mask, 'x'],
                    overlay_points.loc[mask, 'y'],
                    color=colors[celltype],
                    s=30,
                    marker='o',
                    label=celltype
                )
        
        plt.gca().set_aspect(1 / (max(abs(optimally_refined_neighborhood['y'])) / max(optimally_refined_neighborhood['x'])))
        plt.title(f'delta_z_smooth > {round(np.quantile(neighborhoods["delta_z_smooth"], optimal_threshold), 3)}')
        plt.legend(title="Celltype")
        plt.tight_layout()
        plt.savefig(f"{out_dir}/optimal_neighborhood-targets.pdf")
        plt.savefig(f"{out_dir}/optimal_neighborhood-targets.png")
        plt.close()
    
    # Create animated GIFs if requested
    if make_gif:
        try:
            # Create GIF for refinement plots
            refinement_images = []
            for filepath in sorted(glob.glob(f"{out_dir}/refinement/*.png")):
                refinement_images.append(Image.open(filepath))
            
            if refinement_images:
                refinement_images[0].save(
                    f"{out_dir}/refinement.gif",
                    save_all=True,
                    append_images=refinement_images[1:],
                    duration=250,  # 4 fps = 250ms per frame
                    loop=0
                )
            
            # Create GIF for objective function plots
            objective_images = []
            for filepath in sorted(glob.glob(f"{out_dir}/refinement-objective/*.png")):
                objective_images.append(Image.open(filepath))
            
            if objective_images:
                objective_images[0].save(
                    f"{out_dir}/refinement-objective.gif",
                    save_all=True,
                    append_images=objective_images[1:],
                    duration=250,  # 4 fps = 250ms per frame
                    loop=0
                )
        except Exception as e:
            logging.error(f"Error creating GIFs: {e}")
    
    return optimally_refined_neighborhood
