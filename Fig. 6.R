
##########
# mouse
DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .05, ncol = 1, group.by = "celltype_l1")
DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .05, ncol = 1, group.by = "cancer")

input.para = 'Cd74'
input.para = 'H2-Ab1'
input.para = 'H2-Aa'

umap_mat = All_obj_time1@reductions$umap@cell.embeddings
umap_mat = cbind(umap_mat, (All_obj_time1@assays$RNA@data[input.para,]))
colnames(umap_mat)[3] = "para"
umap_mat = data.frame(umap_mat)
quantile(umap_mat$para)
umap_mat_sig = subset(umap_mat, umap_mat$para>quantile(umap_mat$para, 0.99))


umap_mat_large = data.frame(All_obj_time1@reductions$umap@cell.embeddings)

cc <- colorRampPalette(c("white", "#F47E5D", "#CA3D74", "#7F2880", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"

ggplot()+
  geom_point(data = umap_mat_large,  aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.001)+
  geom_point(data =umap_mat_sig, aes(x = UMAP_1, y = UMAP_2, size = para), shape = 1,colour = "black", stroke = 1)+
  geom_point(data =umap_mat_sig, aes(x = UMAP_1, y = UMAP_2, size = para, color = para))+
  labs(x = "UMAP_1", y = "UMAP_2", subtitle = input.para)+
  scale_color_gradientn(colours = cc(100)) +
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  # scale_fill_gradientn(colours = (cc(100))) +
  NULL



library(UCell)
signature1 <- qusage::read.gmt("/www/data/signature/NI_signature_NEU_220504.gmt")
input.mat = All_obj_time1@assays$RNA@data
row.names(input.mat) = toupper(row.names(input.mat) )
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(input.mat, features = signature1, maxRank = 3000, name = "")))
signature.matrix.neu[1:3,1:3]
All_obj_time12 = All_obj_time1
All_obj_time12@meta.data = cbind(All_obj_time12@meta.data, signature.matrix.neu)
FeaturePlot(All_obj_time12, features = "Antigen.Processing.presentation_DIY", pt.size = .1, ncol = 1, 
            cols = c("grey90", "#F47E5D", "#CA3D74", "#7F2880", "#463873"), reduction = "umap", 
            raster = F) & theme(aspect.ratio=1) #& scale_color_viridis_c()





#########3
#mouse volumn
mousemat = read.table("/home/data/data/PC/mouse/mouse2_volumn2_MC38_221208.txt", sep="\t", header=T)
# mousemat = read.table("/home/data/data/PC/mouse/mouse2_volumn2_H1-6_221216.txt", sep="\t", header=T)
# mousemat = read.table("/home/data/data/PC/mouse/mouse2_volumn2_LLC_221222.txt", sep="\t", header=T)
mousemat = subset(mousemat, !is.na(mousemat$Sizea))

mousemat = data.frame(mousemat)

mousemat$volumn = as.numeric(as.character(mousemat$volumn))

mousemat_arranged = mousemat
mousemat_arranged$timepoint = as.Date(mousemat_arranged$timepoint)
mousemat_arranged$implantday = as.Date(mousemat_arranged$implantday)
mousemat_arranged$day = as.numeric(mousemat_arranged$timepoint  - mousemat_arranged$implantday)

mousemat_arranged$group_timepoint = paste(mousemat_arranged$timepoint, mousemat_arranged$group, sep="_")

ggplot(data = mousemat_arranged_select, aes(x = group, y = volumn))+
  ggnewscale::new_scale_color() +
  # geom_bar(stat = "identity", width=0.6, fill = "grey94", color = "black") + #
  geom_boxplot(outlier.shape = NA)+ # 背景色透明化
  geom_point(data = mousemat_arranged_select, aes(fill = group, color = group), size=3, #, shape = group
             position = position_jitterdodge(dodge.width = 0.2), 
             alpha=0.7)+ # 边框线黑色
  theme_bw()+
  theme(aspect.ratio=1)+
  theme(legend.position="none", legend.box = "none", axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 1500, linetype = "longdash") +
  # ylim(0,50)+
  stat_compare_means(comparisons = list(
    c("g1","g2"),c("g1","g7"),c("g7","g2"),
    c("g3","g4"),c("g3","g8"),c("g4","g8"),
    c("g5","g6"),c("g5","g7"),c("g6","g7"),
    c("g6","g8")
  ),method = "t.test",label = "p.signif") + #,label = "p.signif"
  # stat_compare_means(method = "t.test")+ #wilcox.test
  facet_wrap(~cancer + day, ncol = 3, scales = "free") +
  NULL



data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

mousemat_arranged_select2 = data_summary(mousemat_arranged_select, varname="volumn", 
                                         groupnames=c("group", "cancer", "day"))

mousemat_arranged_select2 <- Rmisc::summarySE(mousemat_arranged_select, measurevar="volumn", 
                                              groupvars=c("group", "cancer", "day"))


library(ggpubr)
library(ggplot2)


ggplot(mousemat_arranged_select2, aes(x=day, y=volumn, group=group, color=group)) + 
  geom_line() +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=volumn-se, ymax=volumn+se), width=1,
                position=position_dodge(0.05))+
  theme_classic() +
  # ylim(0,1500) +
  xlim(12,25) +
  scale_color_manual(values=my36colors) +
  NULL




mousemat_arranged_select2_1 = subset(mousemat_arranged_select2, mousemat_arranged_select2$group %in% c("g1","g2","g7"))
p1 = ggplot(mousemat_arranged_select2_1, aes(x=day, y=volumn, group=group, color=group)) + 
  geom_line(size = 1) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=volumn-se, ymax=volumn+se), width=1.5, size = 1,
                position=position_dodge(0.05))+
  theme_classic() +
  ylim(0,800) +
  xlim(11,25) +
  scale_color_manual(values=my36colors[c(1,2,7)]) +
  theme(aspect.ratio=1)+
  NULL
p1


mousemat_arranged_select2_2 = subset(mousemat_arranged_select2, mousemat_arranged_select2$group %in% c("g3","g4","g8"))
p2 = ggplot(mousemat_arranged_select2_2, aes(x=day, y=volumn, group=group, color=group)) + 
  geom_line(size = 1) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=volumn-se, ymax=volumn+se), width=1.5, size = 1,
                position=position_dodge(0.05))+
  theme_classic() +
  ylim(0,1500) +
  xlim(11,25) +
  scale_color_manual(values=my36colors[c(3,4,15)]) + #g8 is color15
  theme(aspect.ratio=1)+
  NULL
p2




mousemat_arranged_select2_3 = subset(mousemat_arranged_select2, mousemat_arranged_select2$group %in% c("g5","g6","g7","g8"))
p3 = ggplot(mousemat_arranged_select2_3, aes(x=day, y=volumn, group=group, color=group)) + 
  geom_line(size = 1) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=volumn-se, ymax=volumn+se), width=1.5, size = 1,
                position=position_dodge(0.05))+
  theme_classic() +
  # ylim(0,1500) +
  xlim(11,25) +
  scale_color_manual(values=my36colors[c(5,6,7,15)]) + #g8 is color11
  theme(aspect.ratio=1)+
  NULL
p3

cowplot::plot_grid(p1,p2,p3, nrow = 1)




