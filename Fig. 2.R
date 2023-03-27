
##########
#too-many-cells
TAN_obj_tmc = All_obj_time1

table(TAN_obj_tmc$cancer)
table(TAN_obj_tmc$type)
TAN_obj_tmc = subset(TAN_obj_tmc, cancer != "Healthy")
TAN_obj_tmc = subset(TAN_obj_tmc, type %in% c("Cancer","Metastasis"))
table(TAN_obj_tmc$cancer)
table(TAN_obj_tmc$type)

input.num = 10000
cellid<-sample(1:ncol(TAN_obj_tmc), input.num, replace=F); length(cellid)
obj_random<-TAN_obj_tmc[,cellid]
dim(obj_random)


obj_random$group_tmc = obj_random$cancer
# obj_random$group_tmc = obj_random$celltype_l1
table(obj_random$group_tmc)



exp4<-data.frame(t(obj_random@reductions$harmony@cell.embeddings))
# exp4<-data.frame(as.matrix(obj_random@assays$RNA@counts[obj_random@assays$RNA@var.features[1:3000],]))
# exp4<-data.frame(as.matrix(obj_random@assays$RNA@scale.data[obj_random@assays$RNA@var.features[1:500],]))
exp4<-tibble::rownames_to_column(exp4, "Gene")
exp4[1:3,1:3]

labels<-data.frame(row.names(obj_random@meta.data), obj_random@meta.data$group_tmc); colnames(labels)<-c("item", "label")
labels[1:3,]
table(labels[,2])

dim(exp4)
dim(labels)

colnames(exp4)<-gsub("-", ".", colnames(exp4))
colnames(exp4)<-gsub("_", ".", colnames(exp4))
colnames(exp4)<-gsub(" ", ".", colnames(exp4))
colnames(exp4)<-gsub("'", ".", colnames(exp4))
labels[,1]<-gsub("_", ".", labels[,1])
labels[,1]<-gsub("-", ".", labels[,1])
labels[,1]<-gsub(" ", ".", labels[,1])
labels[,1]<-gsub("'", ".", labels[,1])

sum(!(colnames(exp4) %in% labels[,1]))

fwrite(exp4, file = "/www/data/PC/treateddata/toomanycells/expr_count.csv", sep = ",")
fwrite(labels, file = "/www/data/PC/treateddata/toomanycells/labels.csv", sep = ",")

#linux 
cd /www/data/PC/treateddata/toomanycells
rm -rf out; mkdir out

#test docker
# sudo docker run -it --rm -v "/www/data/PC:/www/data/PC" gregoryschwartz/too-many-cells:2.0.0.0 -h


docker run -it --rm -v "/www/data/PC:/www/data/PC" gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
--min-size 200 \
--matrix-path /www/data/PC/treateddata/toomanycells/expr_count.csv \
--labels-file /www/data/PC/treateddata/toomanycells/labels.csv \
--draw-collection "PieChart" \
--draw-colors "[\"#5F3D69\", 
\"#E5D2DD\", \"#53A85F\", \"#F1BB72\", \"#F3B1A0\", \"#D6E7A3\", \"#57C3F3\",
\"#E95D5A\", \"#E59DC5\", 
\"#AB3583\", \"#486E88\", \"#23452F\", \"#BD956A\", \"#8C549C\", \"#585658\", 
\"#9FA3A8\", \"#E0D4CA\", \"#C5DEBA\", \"#58A4C3\", \"#E4C755\", \"#F7F398\"]" \
--output /www/data/PC/treateddata/toomanycells/out \
--dendrogram-output ./dendrogram.pdf \
> clusters.csv






#########3
# gene expression


#########
# cytokine
immunomodulator = read.table("/home/wuyingcheng/data/reference/immunomodulator.txt", header=F, sep="\t")
table(immunomodulator$V2)
immunomodulator = subset(immunomodulator, immunomodulator$V2 != "MHC")

#MHC
immunomodulator = subset(immunomodulator, immunomodulator$V2 == "MHC")


# ISG
signature1 <- qusage::read.gmt("/www/data/signature/c5.go.v7.4.symbols.gmt")
names(signature1)[names(signature1) %like% "INTERFERON"]
isggene = c(signature1[["GOBP_INTERFERON_ALPHA_PRODUCTION"]],
            signature1[["GOBP_INTERFERON_BETA_PRODUCTION"]],
            signature1[["GOBP_INTERFERON_GAMMA_PRODUCTION"]])

DimPlot(All_obj_time1, group.by = c("celltype_l1"), ncol = 1, label = T)

unique((All_obj_time1$celltype_l1))
input.order = c("HLA-DR+CD74+","IL1R2+NFKBIA+","CXCL8+","IFIT1+ISG15+","CSF3R+","S100A12+","TXNIP+","MMP9+CEACAM8+","VEGFA+SPP1+")
Idents(All_obj_time1) = All_obj_time1$celltype_l1
Idents(All_obj_time1) <- factor(Idents(All_obj_time1), levels= input.order)
All_obj_time1$celltype_l1 <- factor(All_obj_time1$celltype_l1, levels= input.order)
unique(Idents(All_obj_time1))

library(COSG)
marker_cosg <- cosg(All_obj_time1, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
marker_cosg_wide = marker_cosg$names
all.gene = as.character(as.matrix(marker_cosg_wide))
all.gene.select = all.gene[all.gene %in% immunomodulator$V3]

library(scRNAtoolVis)

Idents(All_obj_time1) = All_obj_time1$celltype_l1
DotPlot(All_obj_time1, features = unique(all.gene.select))+ ggplot2::coord_flip() +RotatedAxis()


gene_cell_exp <- AverageExpression(All_obj_time1, features = unique(all.gene.select),slot = 'data', group.by = "celltype_l1") 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

pheatmap(gene_cell_exp,
         show_rownames = T, show_colnames = T,
         #color = pal,
         cluster_rows = F, cluster_cols = F,
         #clustering_method = input.clustering , #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         scale = "row",
         #annotation_col = annotation_row,
         #breaks = breaks,
         #annotation_row = annotation_row,
         #clustering_distance_cols = input.clusterdistance,
         #annotation_row = annotation_row,
         #cellwidth = 3, cellheight = 3,
         # color = (cc(100)),
         border_color = NA
)




gene_cell_exp_tmp = gene_cell_exp
gene_cell_exp_tmp$rank = row.names(gene_cell_exp_tmp)
gene_cell_exp_long = reshape2::melt(gene_cell_exp_tmp, id.vars = "rank")

gene_cell_exp_long = DotPlot(All_obj_time1, features = (unique(all.gene.select)))$data


input.order = rev(c("HLA-DR+CD74+","IL1R2+NFKBIA+","CXCL8+","IFIT1+ISG15+","CSF3R+","S100A12+","TXNIP+","MMP9+CEACAM8+","VEGFA+SPP1+"))
gene_cell_exp_long$id = as.character(gene_cell_exp_long$id)
gene_cell_exp_long$id <- factor( gene_cell_exp_long$id, levels= input.order)

ggplot(gene_cell_exp_long, aes(features.plot, id)) + 
  geom_point(aes(size = pct.exp, fill =  avg.exp.scaled, color = avg.exp.scaled), shape = 21, colour = "black")+
  scale_fill_gradientn(colours = (cc(100))) +
  NULL






