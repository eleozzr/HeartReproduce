#infer pseudotime for fibroblast from DCM heart
library(monocle3)
library(Seurat)
library(grid)
library(gglot2)
library(dplyr)

obj=readRDS("../obj_fb_DCM_processed_label.rds")
obj=NormalizeData(obj)
pData<-obj@meta.data
fData <- data.frame(gene_short_name = rownames(obj), row.names = rownames(obj))
#Step1 Creat monocle3 object
cds <- new_cell_data_set(obj@assays$originalexp@counts,
                         cell_metadata =  pData,
                         gene_metadata =  fData)
data<-preprocess_cds(cds,num_dim=50,method="PCA")
#plot_pc_variance_explained(data)
data_align <- align_cds(data, alignment_group = "sample",alignment_k=20)
data_align<-reduce_dimension(data_align,umap.min_dist=0.20,umap.n_neighbors = 15,umap.fast_sgd=F,reduction_method="UMAP",preprocess_method="Aligned")

#Cluster using different resolution
data_align_cluster <- cluster_cells(data_align,resolution=0.0003,random_seed = 2307)
data_align_cluster$clusters1 <- data_align_cluster@clusters$UMAP$clusters
data_align_cluster <- monocle3::cluster_cells(data_align_cluster,resolution=5e-4,random_seed=2307)
data_align_cluster$clusters2 <- data_align_cluster@clusters$UMAP$clusters
data_align_cluster <- monocle3::cluster_cells(data_align_cluster,resolution=8e-4,random_seed=2307)
data_align_cluster$clusters3 <- data_align_cluster@clusters$UMAP$clusters
data_align_cluster <- learn_graph(data_align_cluster,close_loop = F,learn_graph_control=list(minimal_branch_len=15,geodesic_distance_ratio=1/2,rann.k=50))
#plot_cells(data_align_cluster, color_cells_by = "cluster",cell_size=1)
#plot_cells(data_align_cluster,cell_size=1,genes="PKNOX2",label_cell_groups=FALSE)
#data_align_cluster<- order_cells(data_align_cluster,verbose=F,root_cells=colnames(data_align_cluster)[data_align_cluster$subfb=="MSC"])

#choose root cells
mst <- principal_graph(data_align_cluster)$UMAP
cell_ids <- colnames(data_align_cluster)[data_align_cluster$subfb ==  "MSC"]
closest_vertex <- data_align_cluster@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(data_align_cluster), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(data_align_cluster)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
#inference psedutotime 
data_align_cluster <- order_cells(data_align_cluster, root_pr_nodes = root_pr_nodes)
#plot_cells(cds, color_cells_by = "pseudotime",cell_size=1)
#save object

save(data_align_cluster,mst,file="cds_fb_DCM_monocle3_processed.RData")

logcounts(data_align_cluster)=counts(data_align_cluster)
obj_seurat=as.Seurat(data_align_cluster)
saveRDS(obj_seurat,file="obj_fb_DCM_monocle3_processed.rds")

