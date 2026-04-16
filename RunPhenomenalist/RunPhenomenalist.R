RunPhenomenalist=function(segmentation_file,failed.markers=NULL,nuclear.markers=NULL,HALO=T,mask.only=NULL,out_dir=getwd(),clustering_res=seq(5,7),classifier_label=NULL,min.cells=10,else_cytoplasm=F,max.cells=1e5,phenotyping_template=NULL){
  suppressPackageStartupMessages({
    library(phenomenalist)
    library(tidyverse)
    library(glue)
    library(cowplot)
    library(ggsci)
  })

  # Locate phenomenalist-utils.R (local .mod function variants + prepare_mask_inputs
  # helpers that are not part of the public phenomenalist package). Search order:
  #   1. PHENOMENALIST_DIR env var (set by run-phenomenalist.R)
  #   2. Same directory as the currently running script
  #   3. Current working directory
  # Skip if already sourced.
  if (!exists("plot_heatmap.mod_v1", mode = "function")) {
    .pheno_dir <- Sys.getenv("PHENOMENALIST_DIR", unset = "")
    .candidates <- c(
      if (nzchar(.pheno_dir)) file.path(.pheno_dir, "phenomenalist-utils.R"),
      file.path(getwd(), "phenomenalist-utils.R")
    )
    .utils <- Find(file.exists, .candidates)
    if (is.null(.utils)) {
      stop("phenomenalist-utils.R not found. Set PHENOMENALIST_DIR or run from the ",
           "RunPhenomenalist/ directory.")
    }
    source(.utils)
  }

  # Default cols to skip for HALO:
  basic_skip_cols = "Blank|blank|DAPI"
  
  data_loc=segmentation_file
  
  ## Default cols to skip for HALO:
  # here if using
  if(is.null(nuclear.markers)){
    halo_skip_cols = "Classification|Completeness|Cytoplasm|Nucleus|Membrane|DAPI"
    if(length(failed.markers)>0){
      halo_skip_cols = "Classification|Completeness|Cytoplasm|Nucleus|Membrane|DAPI"
      halo_skip_cols=paste0(halo_skip_cols,paste(failed.markers,collapse = '|'),collapse = '|')
      
    }
  }else{
    halo_skip_cols = "Classification|Completeness|Membrane|"
    if(length(failed.markers)>0){
      halo_skip_cols=paste0(halo_skip_cols,paste(failed.markers,collapse = '|'),collapse = '|')
      
    }
  }
  
  
  
  
  
  
  
  
  labels.files=unlist(strsplit(data_loc,'[/]'))[length(unlist(strsplit(data_loc,'[/]')))]
  label=unlist(strsplit(labels.files,'.csv'))[1]
  
  
  
  # by default NULL:
  expression.columns=phenomenalist.preprocess(segmentation_file,failed.markers = failed.markers,nuclear.markers = nuclear.markers,HALO = HALO,else.cytoplasm = else_cytoplasm)
  message(expression.columns)
  print(expression.columns)
  if(length(expression.columns) == 0){
	
	expression.columns=NULL

  }  
  if(!HALO){
    skip_cols=paste0(c(basic_skip_cols,'Ch'),collapse = '|')
    if(length(failed.markers)>0){
      skip_cols=paste0(skip_cols,'|')
      skip_cols=paste0(skip_cols,paste(failed.markers,collapse = '|'),collapse = '|')
    }
    
    # if(length(data_loc)==1){
    #   data_csv=data_loc
    #   print(data_csv)
    #   x <- data.table::fread(data_csv, stringsAsFactors = FALSE, data.table = FALSE)
    #   expression.columns=colnames(x)[seq(6,dim(x)[2])]
    #   expression.columns=expression.columns[!str_detect(expression.columns,skip_cols)]
    # }
  }else{
    skip_cols=halo_skip_cols
  }
  
  
  # clustering resolutions:
  resolutions = clustering_res
  # which subset to store (by default all)
  idx=seq(length(resolutions))
  
  print(resolutions)
  
  # Where output dirs should be housed:
  workDir=out_dir
  
  if (!dir.exists(workDir)) { dir.create(workDir) }
  setwd(workDir)
  
  out.label=NULL
  
  data_csv=segmentation_file
  message(glue('Reading file: {data_csv}'))
  
  
  message(glue('Skipping columns: {skip_cols}'))
  
  out_dir = "./out-phenomenalist"
  if (!dir.exists(out_dir)) { dir.create(out_dir) }
  
  out_dir = paste0(out_dir, "/", label)
  message(out_dir)
  
  if('lock.rds' %in% list.files(out_dir)){
    spe = readRDS(glue('{out_dir}/spe.rds'))
    spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
    names(colData(spe))
    saveRDS(spe,glue(out_dir,'/spe.rds'))
    
    cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
    
    for(i in seq(length(cluster_cols))){
      
      plot_heatmap.mod_v1(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
    }
    
    res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
    prepare_mask_inputs(spe=spe,out_dir = out_dir,res = max(as.numeric(res_found[,19])),mask.only = mask.only,label=label)
    
    
  }else{
    spe <- create_object.mod(data_csv, skip_cols = skip_cols, transformation = "z", out_dir = out_dir,expression_cols = expression.columns,classifier.label = classifier_label,min.cells = min.cells,max.cells = max.cells)
    
    if ("tile_num" %in% names(colData(spe))) plot_spatial(spe, color_by = "tile_num",out_dir = out_dir)
    
    if ("Classifier_Label" %in% names(colData(spe))) plot_spatial(spe, color_by = "Classifier_Label",out_dir = out_dir)
    
    spatial_dir = paste0(out_dir, "/spatial-expression")
    dir.create(spatial_dir, showWarnings = FALSE)
    plot_spatial(spe, color_by = names(spe), out_dir = spatial_dir)
    plot_spatial(spe, color_by = names(spe), smooth = TRUE, out_dir = spatial_dir)
    
    
    umap_dir = paste0(out_dir, "/UMAP-expression")
    dir.create(umap_dir, showWarnings = FALSE)
    plot_dr(spe, dr = "UMAP", color_by = names(spe), out_dir = umap_dir)
    plot_dr(spe, dr = "UMAP", color_by = names(spe), smooth = TRUE, out_dir = umap_dir)
    
    phenotyping_dir=paste0(out_dir,"/phenotyping_template_results")
    if(!is.null(phenotyping_template)){
	
	template=read.csv(phenotyping_template)
	dir.create(phenotyping_dir)
    }
    
    spe = cluster.mod(spe, resolution = resolutions, out_dir = out_dir)
    print(names(colData(spe)))
    saveRDS(spe,glue(out_dir,'/spe.rds'))
    
    cluster_cols = stringr::str_subset(names(colData(spe)), "cluster_leiden")
    
    for(i in seq(length(cluster_cols))){
      
      plot_heatmap.mod_v1(spe, group_by = cluster_cols[i],auto = T,out_dir = out_dir,segment=T)
      if(!is.null(phenotyping_template)){
	spe = assign_celltype_with_template(obj=spe,phenotyping_template=template,cluster=cluster_cols[i])
	plot_heatmap.mod_v1(spe,group_by=glue('{cluster_cols[i]}_annotations_template'),out_dir=phenotyping_dir)
	plot_dr(spe, dr = "UMAP", color_by = glue('{cluster_cols[i]}_annotations_template'), smooth = TRUE, out_dir = phenotyping_dir)
	plot_spatial(spe, color_by = glue('{cluster_cols[i]}_annotations_template'), smooth = TRUE, out_dir = phenotyping_dir)
	spatial_obj=data.frame(spatialCoords(spe))
	spatial_obj$cluster=obj[[glue('{cluster_cols[i]}_annotations_template')]]
	write.csv(spatial_obj,glue('{phenotyping_dir}/{cluster_cols[i]}_annotations_template.csv'))
      }

    }
    
    res_found=do.call('rbind',strsplit(cluster_cols,'[cluster_leiden_res]'))
    for(res in seq(dim(res_found)[1])){
   	 prepare_mask_inputs(spe=spe,out_dir = out_dir,res =as.numeric(res_found[res,19]),mask.only = mask.only,label=label)
    }
  }
  
  
  
  
}
