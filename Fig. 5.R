
##########
# correlation between cell types (matched RNA-seq)
PC_immunemat = data.frame(cbind(t(PC_xCell), PC_NEU))
selectcell = c(colnames(PC_NEU),"B.cells","CD4..memory.T.cells","CD4..naive.T.cells","CD4..T.cells","CD4..Tcm","CD4..Tem","CD8..naive.T.cells","CD8..T.cells","CD8..Tcm","CD8..Tem","Class.switched.memory.B.cells","DC","Endothelial.cells","Epithelial.cells","Fibroblasts","Macrophages","Macrophages.M1","Macrophages.M2","Memory.B.cells","Monocytes", "naive.B.cells","NK.cells","NKT","pDC","Plasma.cells","pro.B.cells","Tgd.cells","Th1.cells","Th2.cells","Tregs")
PC_immunemat = PC_immunemat[,selectcell]

cor_mat = cor(PC_immunemat, method = "spearman")
# cor_mat[cor_mat < 0] = 0

library(pheatmap)

cc <- colorRampPalette(c("white", "#fc58a6"))

pheatmap(cor_mat, show_rownames = T, show_colnames = T,
         cluster_rows = T, cluster_cols = T,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         scale = "none",
         cellwidth = 7, cellheight = 7,
         color = (cc(100)),
         border_color = NA
)








#########
#HLA
mtx.mat_l1 = data.frame(mtx.mat_l1)
mtx.mat_l1 = subset(mtx.mat_l1, rowMeans(mtx.mat_l1) > 0)
mtx.mat_l2 = data.frame(mtx.mat_l2)
mtx.mat_l2 = subset(mtx.mat_l2, rowMeans(mtx.mat_l2) > 0)

All_obj_time1_sub@meta.data = cbind(All_obj_time1_sub@meta.data, t(mtx.mat_l1))
All_obj_time1_sub@meta.data = cbind(All_obj_time1_sub@meta.data, t(mtx.mat_l2))

input.features_l1 = row.names(mtx.mat_l1)
input.features_l2 = row.names(mtx.mat_l2)
FeaturePlot(All_obj_time1_sub, features = input.features_l1, pt.size = 0.1, ncol = 5, reduction = "umap", raster = F)
Nebulosa::plot_density(All_obj_time1_sub, input.features_l1, size = 1)


stat_hla = data.frame(type = row.names(mtx.mat_l1), 
                      mean = rowMeans(mtx.mat_l1),
                      num = rowSums(mtx.mat_l1 > 0)
)
stat_hla = stat_hla[order(stat_hla$num, decreasing = T),]
stat_hla[1:10,]
stat_hla$lognum = log10(stat_hla$num)

stat_hla = subset(stat_hla, stat_hla$type %like% ":")
stat_hla_1 = subset(stat_hla, !stat_hla$type %like% "HLA-D")
stat_hla_2 = subset(stat_hla, stat_hla$type %like% "HLA-D")


stat_hla_1_top = stat_hla_1[1:10,]
stat_hla_2_top = stat_hla_2[1:10,]

ggplot(stat_hla_1, aes(y = lognum, x = mean)) + 
  geom_point(aes(size = lognum), color = "black",shape= 21, stroke = 1) +
  geom_point(aes(size = lognum, color = mean), shape= 16) +
  ggrepel::geom_text_repel(stat_hla_1_top, mapping = aes(label = type), show.legend = F,max.overlaps =10000,
                           segment.size = 0.2, segment.alpha = 0.5) +
  # scale_color_manual(values=my36colors) + #my36colors
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradientn(colours = cc(100)) +
  ylim(0.9,4.2)+
  theme_bw()+
  theme(aspect.ratio=.5,axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "right")+
  NULL


ggplot(stat_hla_2, aes(y = lognum, x = mean)) + 
  geom_point(aes(size = lognum), color = "black",shape= 21, stroke = 1) +
  geom_point(aes(size = lognum, color = mean), shape= 16) +
  ggrepel::geom_text_repel(stat_hla_2_top, mapping = aes(label = type), show.legend = F,max.overlaps =10000,
                           segment.size = 0.2, segment.alpha = 0.5) +
  # scale_color_manual(values=my36colors) + #my36colors
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradientn(colours = cc(100)) +
  ylim(0.8,4)+
  theme_bw()+
  theme(aspect.ratio=.5,axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "right")+
  NULL


