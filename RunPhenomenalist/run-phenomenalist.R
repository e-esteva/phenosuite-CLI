source('/gpfs/data/abl/tric/segmentation/CODEX-Pipeline/RunPhenomenalist/RunPhenomenalist.R')

library(glue)
library(MatrixGenerics)
library(SpatialExperiment)

arguments=commandArgs(TRUE)
# 1: segmentation_file
# 2: failed.markers
# 3: nuclear.markers
# 4: HALO
# 5: out_dir
# 6: clustering_res
# 7: classifier_label
# 8: max cells
# 9: phenotyping template

arguments=data.frame(t(arguments))

names(arguments)=c('segmentation_file','failed.markers','nuclear.markers','HALO','out_dir','clustering_res','classifier_label','max.cells','phenotyping_template')
print(arguments)
print('initiating')

print(arguments$max.cells)
print(length(arguments))

phenotyping_template=arguments$phenotyping_template

HALO=as.logical(arguments$HALO)
message(glue('HALO: {HALO}'))

if(arguments$out_dir == "0"){
        arguments$out_dir=getwd()
}
if(arguments$clustering_res=="0"){
        clustering_res=seq(5,7)

}else{
	res=unlist(strsplit(as.character(arguments$clustering_res),'[,]'))
        if(length(res) == 2){
		clustering_res=seq(as.numeric(res[1]),as.numeric(res[2]))
	}else{
		clustering_res=as.numeric(res)
	}
}
print('clustering resolution: ')
print(clustering_res)


if(arguments$failed.markers=='0'){
        failed.markers=NULL

}else{
      	failed.markers=unlist(strsplit(as.character(arguments$failed.markers),'[,]'))
}


if(arguments$nuclear.markers=='0'){
        nuclear.markers=NULL

}else{
      	nuclear.markers=unlist(strsplit(as.character(arguments$nuclear.markers),'[,]'))
}


if(arguments$classifier_label=='0' || arguments$classifier_label=='NULL'){
        classifier_label=NULL

}else{
      	classifier_label=unlist(strsplit(as.character(arguments$classifier_label),'[,]'))
}

if(arguments$phenotyping_template==''){
	phenotyping_template=NULL

}
arguments$max.cells=as.numeric(arguments$max.cells)

lapply(list.files('/gpfs/data/abl/tric/segmentation/CODEX-Pipeline/phenomenalist/R/'),function(x) source(glue('/gpfs/data/abl/tric/segmentation/CODEX-Pipeline/phenomenalist/R/{x}')))

RunPhenomenalist(segmentation_file=arguments$segmentation_file,failed.markers=failed.markers,nuclear.markers=nuclear.markers,HALO=HALO,out_dir=arguments$out_dir,clustering_res=clustering_res,classifier_label=classifier_label,max.cells=arguments$max.cells,min.cells=10,phenotyping_template=phenotyping_template)


