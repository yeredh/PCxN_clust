rm(list=ls())
options(stringsAsFactors = F)

# load function to transform correlations into distance matrix
source("src/pcxn_clust_functions.R")

# ==== PCxN: CP ====
# Data frame with pathway correlations for the Canonical Pathways
pcxn = readRDS("data/pathCor_CPv5.1_dframe.RDS")


# Default filters:
# p-value < 0.05 (significant correlation coefficients)
pcxn = pcxn[pcxn$p.Adjust < 0.05,]
# |PathCor| > 0.05 (pathway correlation magnitude)
pcxn = pcxn[abs(pcxn$PathCor) > 0.05,]



# ==== KEGG AD Top 10 neighbors ====
# get top 10 neighbors of KEGG_ALZHEIMERS_DISEASE
pcxn_ad=getTopNeighbors(path_name="KEGG_ALZHEIMERS_DISEASE",pcxn=pcxn,top_n=10)




library(ape)
library(rafalib)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
# ==== Clustering: Correlation Matrix ====
# get correlation matrix
pcxn_cor = pcxn2cormat(dat=pcxn_ad)
# Euclidean distance, complete method
hc_cor = hclust(d=dist(x=pcxn_cor,method="euclidean"), method="complete")

# Heatmap
hm_cor = Heatmap(
  matrix=pcxn_cor,cluster_rows = hc_cor,cluster_columns = hc_cor,name="PathCor",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(9, "cm"),
  column_names_max_height = unit(9, "cm"),
  col = colorRamp2(c(min(pcxn_cor[lower.tri(pcxn_cor)]), 0, max(pcxn_cor[lower.tri(pcxn_cor)])), c("blue", "white", "red"))
)
draw(hm_cor,heatmap_legend_side = "left")

# Dendogram
plot(as.phylo(hc_cor),cex=0.75,font=2)


# ==== Clustering: d=1-|PathCor| ====
# complete, d= 1-|PathCor|
hc_out = hclust(d=pcxn2dist(dat=pcxn_ad),method="complete")

# Heatmap
hm_dist = Heatmap(
  matrix=pcxn_cor,cluster_rows = hc_out,cluster_columns = hc_out,name="PathCor",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(9, "cm"),
  column_names_max_height = unit(9, "cm"),
  col = colorRamp2(c(min(pcxn_cor[lower.tri(pcxn_cor)]), 0, max(pcxn_cor[lower.tri(pcxn_cor)])), c("blue", "white", "red"))
)
draw(hm_dist,heatmap_legend_side = "left")

# Dendogram
plot(as.phylo(hc_out),cex=0.75,font=2)





library(RColorBrewer)
palette(brewer.pal(n=8, name="Dark2"))
# ==== Plots using d=1-|PathCor| ==== 
# cut tree to define 4 clusters
clustMember = cutree(hc_out, k=4)

# Heatmap: assign color to clusters
hm_dist = Heatmap(
  matrix=pcxn_cor,
  cluster_rows = hc_out,
  cluster_columns = hc_out,name="PathCor",
  row_names_gp = gpar(fontsize = 7,col=clustMember,font=2),
  column_names_gp = gpar(fontsize = 7,col=clustMember,font=2),
  row_names_max_width = unit(9, "cm"),
  column_names_max_height = unit(9, "cm"),
  col = colorRamp2(c(min(pcxn_cor[lower.tri(pcxn_cor)]), 0, max(pcxn_cor[lower.tri(pcxn_cor)])), c("blue", "white", "red"))
)
draw(hm_dist,heatmap_legend_side = "left")

# Dendogram: assign color to clusters
plot(as.phylo(hc_out),cex=0.75,tip.col = clustMember,font=2)

library(igraph)
library(plotrix)
# ==== Plot: Graph =====
# plot network based on pathway correlations
# - node color: matches the clusters from the hclust output
# - edge color: green for PathCor >= 0, red for PathCor < 0
# - edge width: proportional to |PathCor|
# - layout: based on distance used in hierarchical clustering

# estimate distance
pcxn_ad$D = 1- abs(pcxn_ad$PathCor)
# create graph
g = graph_from_data_frame(d=pcxn_ad[,c("Pathway.A","Pathway.B","PathCor","D")],directed = F)


# set up color palette for clusters
col_pal = clustMember[names(V(g))]
col_pal[!grepl("[a-zA-Z]+",col_pal)] = palette()[as.numeric(col_pal[!grepl("[a-zA-Z]+",col_pal)])]

# plot graph
set.seed(6)
my_lay = layout_with_fr(graph=g,weights=E(g)$D)
plot(
  g,
  layout=my_lay,
  rescale=T,
  vertex.color=col_pal,
  vertex.label=gsub("_"," \n",names(V(g))),
  vertex.label.cex=0.55,
  vertex.label.font=2,
  vertex.label.color="black",
  vertex.frame.color=NA,
  edge.color=ifelse(E(g)$PathCor >= 0,"forestgreen","red"),
  edge.width=rescale(abs(E(g)$PathCor),c(0,1))*2
)
