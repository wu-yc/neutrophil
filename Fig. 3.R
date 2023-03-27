
##########
#pseudotime and Ro/e

load(file = "/www/data/PC/treateddata/All_obj_time1.rda")



meta.tb = All_obj_time1@meta.data
meta.tb$meta.cluster = meta.tb$celltype_l1
meta.tb$loc = meta.tb$type

table(meta.tb$loc)
meta.tb$loc[meta.tb$loc == "Metastasis"] = "Cancer"

OR.NEU.list <- do.tissueDist(cellInfo.tb = meta.tb,
                             out.prefix=sprintf("test",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

OR.NEU.list_p = OR.NEU.list$p.dist.tb
OR.NEU.list_heat = OR.NEU.list$OR.dist.mtx
OR.NEU.list_bin = OR.NEU.list$freq.dist.bin

OR.NEU.list_heat_rmoutlier= OR.NEU.list_heat


orderedcell = pheatmap::pheatmap(OR.NEU.list_heat_rmoutlier,
                                 show_colnames = T,
                                 scale = "none",
                                 cluster_rows = T, 
                                 cluster_cols = T, 
                                 color = (cc(100)),
                                 border_color = NA,
                                 cellwidth = 10, cellheight = 10,
                                 clustering_distance_cols = "correlation",
)



all.row = rownames(OR.NEU.list_heat_rmoutlier[orderedcell$tree_row[["order"]],])
all.col = colnames(OR.NEU.list_heat_rmoutlier[,orderedcell$tree_col[["order"]]])

OR.NEU.list_heat_rmoutlier2 = (OR.NEU.list_heat_rmoutlier)
roelong = reshape2::melt(OR.NEU.list_heat_rmoutlier2)
colnames(roelong) = c("name","variable","value")

roelong$value = as.numeric(as.character(roelong$value))
roelong$value = round(roelong$value, 2)
roelong$value[roelong$value ==0] = 0


roelong$name <- as.character(roelong$name)
roelong$name <- factor(roelong$name, levels= rev(all.row))

roelong$variable <- as.character(roelong$variable)
roelong$variable <- factor(roelong$variable, levels= all.col)

cc = colorRampPalette((brewer.pal(n = 7,  name = "OrRd")))	
ggplot(roelong, aes(variable, name)) + 
  geom_tile(aes(fill = value, fill = value),size=1)+
  # scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  # scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = my12colors[1])+ # "#4195C4"# ("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
  scale_fill_gradient2(low = "white", high = cc(100)[60])+ # "#4195C4"# ("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
  geom_text(aes(label=value),col ="black",size = 3)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legen
  labs(fill =paste0("RED: NEU signature",""))


roelong_cancer = subset(roelong, roelong$variable == "Cancer")


meta_tmp = All_obj_time1@meta.data
row.names(meta_tmp) = gsub("-", "_", row.names(meta_tmp))
row.names(meta_tmp) = gsub("\\.", "_", row.names(meta_tmp))


sctour = read.csv("/www/data/PC/treateddata/All_obj_time1_sctour_ptime.csv")
sctour$X = gsub("-", "_", sctour$X)
sctour$X = gsub("\\.", "_", sctour$X)

sum(row.names(meta_tmp) %in% sctour$X)
dim(sctour)

meta_tmp$X = row.names(meta_tmp)
meta_tmp = left_join(meta_tmp, sctour, by = "X")



vlndata = meta_tmp
vlndata = subset(vlndata, !is.na(vlndata$ptime))
ggplot(vlndata , aes(x=reorder(celltype_l1, ptime, FUN=median, y = ptime), y = ptime,fill = celltype_l1)) + 
  geom_boxplot(outlier.shape=NA)+ 
  scale_fill_manual(values=my12colors) +
  #geom_smooth(aes(group = 1))+
  theme_bw()+theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1)) +
  NULL


#get median pseudotime
all.celltype_l1 = unique(meta_tmp$celltype_l1)
median_pseudo = c()
for (j in 1:length(all.celltype_l1)){
  meta_tmp2 = subset(meta_tmp, celltype_l1 == all.celltype_l1[j])
  median_pseudo = rbind(median_pseudo, c(median(meta_tmp2$ptime, na.rm = T), all.celltype_l1[j]))
  
}
median_pseudo = data.frame(median_pseudo)
colnames(median_pseudo) = c("ptime","celltype_l1")


roelong_cancer2 = roelong_cancer[,c("name","value")]
colnames(roelong_cancer2) = c("celltype_l1", "OR")

roeptime_df = left_join(median_pseudo, roelong_cancer2, by = "celltype_l1")



roeptime_df$ptime = as.numeric(as.character(roeptime_df$ptime))

ggplot(roeptime_df, aes(x = ptime, y = OR)) + 
  geom_point(aes(size = OR, color = celltype_l1)) +
  theme_bw()+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "right")+
  NULL











