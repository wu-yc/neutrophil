
##########
#consensus abundance

neu_median_mat = c()
all.cancer = as.character(all.cancer)
for (j in 1:length(all.cancer)){
  neu_comb_sub = subset(neu_comb, neu_comb$project_id == all.cancer[j])
  neu_median_mat = rbind(neu_median_mat,
                         c(mean(neu_comb_sub$Neutrophils_xCell, na.rm = T),
                           mean(neu_comb_sub$Neutrophils_MCPcounter, na.rm = T),
                           mean(neu_comb_sub$Neutrophils_quantiseq, na.rm = T),
                           mean(neu_comb_sub$Neutrophils_Immunity, na.rm = T),
                           all.cancer[j])
  )
}
neu_median_mat = data.frame(neu_median_mat)
colnames(neu_median_mat) = c("xCell", "MCPcounter", "quantiseq", "Immunity", "cancer")
row.names(neu_median_mat) = neu_median_mat$cancer
neu_median_mat$cancer = NULL

neu_median_mat[,1] = as.numeric(as.character(neu_median_mat[,1]))
neu_median_mat[,2] = as.numeric(as.character(neu_median_mat[,2]))
neu_median_mat[,3] = as.numeric(as.character(neu_median_mat[,3]))
neu_median_mat[,4] = as.numeric(as.character(neu_median_mat[,4]))



rank(neu_median_mat[,1])


neu_median_mat$xCellrank = rank(neu_median_mat$xCell)
neu_median_mat$MCPcounterrank = rank(neu_median_mat$MCPcounter)
neu_median_mat$quantiseqrank = rank(neu_median_mat$quantiseq)
neu_median_mat$Immunityrank = rank(neu_median_mat$Immunity)

neu_median_mat$consensusrank = rowMeans(neu_median_mat[,c("xCellrank","MCPcounterrank","quantiseqrank")])
neu_median_mat <- neu_median_mat[order(neu_median_mat$consensusrank, decreasing = T),]

neu_median_mat_heat = t(neu_median_mat[,c("consensusrank","MCPcounter","quantiseq", "xCell")])
# neu_median_mat_heat = t(neu_median_mat[,c("consensusrank","MCPcounterrank","quantiseqrank","xCellrank")])



cc <- colorRampPalette(c("#352a86", "#095cd8", "#46b896", "#e7ba4a", "#f8fa0d"))
# cc <- colorRampPalette(c("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
pheatmap::pheatmap(neu_median_mat_heat,
                   show_colnames = T,
                   scale = "row",
                   cluster_rows = F, 
                   cluster_cols = F, 
                   color = (cc(100)),
                   border_color = NA,
                   clustering_distance_cols = "correlation",
                   # clustering_method = "centroid", #'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                   # width = 6,
                   # height = 2.5,
                   # annotation_col = annotation_row
                   # filename = file.path(savePath, paste0("NMF-", rank, "-heatmap.png")),
                   # treeheight_row = 10,
                   # treeheight_col = 14
)


load(file = "/www/data/reference/TCGA/tcga_xcell.rda")
load(file = "/www/data/reference/TCGA/tcga_mcp.rda")
load(file = "/www/data/reference/TCGA/tcga_quantiseq.rda")

tcga_xcell$tmp = substr(row.names(tcga_xcell), 1, 15); tcga_xcell<-tcga_xcell[!duplicated(tcga_xcell$tmp), ]; row.names(tcga_xcell) = tcga_xcell$tmp; tcga_xcell$tmp = NULL
tcga_mcp$tmp = substr(row.names(tcga_mcp), 1, 15); tcga_mcp<-tcga_mcp[!duplicated(tcga_mcp$tmp), ]; row.names(tcga_mcp) = tcga_mcp$tmp; tcga_mcp$tmp = NULL
tcga_quantiseq$tmp = substr(row.names(tcga_quantiseq), 1, 15); tcga_quantiseq<-tcga_quantiseq[!duplicated(tcga_quantiseq$tmp), ]; row.names(tcga_quantiseq) = tcga_quantiseq$tmp; tcga_quantiseq$tmp = NULL

load(file = "/www/data/reference/TCGA/tcga_obj.rda")

sample_intersect = intersect(colnames(tcga_obj), row.names(tcga_xcell))
sample_intersect = intersect(sample_intersect, row.names(tcga_mcp))
sample_intersect = intersect(sample_intersect, row.names(tcga_quantiseq))

tcga_obj_sub = subset(tcga_obj, cells = sample_intersect)
sample_intersect = colnames(tcga_obj_sub)

tcga_xcell = tcga_xcell[sample_intersect,]
tcga_mcp = tcga_mcp[sample_intersect,]
tcga_quantiseq = tcga_quantiseq[sample_intersect,]

tcga_obj_sub$NEU_xcell = tcga_xcell$Neutrophils
tcga_obj_sub$NEU_mcp = tcga_mcp$Neutrophils
tcga_obj_sub$NEU_quantiseq = tcga_quantiseq$Neutrophils


colnames(tcga_obj_sub@meta.data)
tcga_obj_sub$project_id[tcga_obj_sub$project_id == ""] = NA

tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, !is.na(tcga_obj_sub@meta.data$project_id))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, (tcga_obj_sub@meta.data$project_id != "TCGA-LAML"))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, (tcga_obj_sub@meta.data$project_id != "TCGA-DLBC"))))
tcga_obj_sub = subset(tcga_obj_sub, cells = row.names(subset(tcga_obj_sub@meta.data, !(tcga_obj_sub@meta.data$project_id %like% "TARGET"))))



tcga_obj_sub = SCTransform(tcga_obj_sub)
tcga_obj_sub <- RunUMAP(tcga_obj_sub, return.model = TRUE, dims = 1:50)
tcga_obj_sub <- RunTSNE(tcga_obj_sub, return.model = TRUE, dims = 1:50)


DimPlot(tcga_obj_sub, reduction = "tsne", label = T, pt.size = .1, 
        group.by = "project_id",  raster = F)
p1 = DimPlot(tcga_obj_sub, reduction = "tsne", label = T, pt.size = .1, 
             group.by = "project_id",  raster = F)

FeaturePlot(tcga_obj_sub, features = c("NEU_xcell","NEU_mcp","NEU_quantiseq"), pt.size = .01, ncol = 3, reduction = "tsne") & theme_bw() & theme(aspect.ratio=1)

cc <- colorRampPalette(c("#352a86", "#095cd8", "#46b896", "#e7ba4a", "#f8fa0d"))
# cc <- colorRampPalette(c("#000075", "#9c0ef0", "#fc58a6", "#ffa955", "#ffff60")) #matched 2
p2 = Nebulosa::plot_density(tcga_obj_sub, "NEU_xcell","NEU_mcp","NEU_quantiseq", size = 0.5, reduction = "tsne") & theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank()) &
  scale_color_gradientn(colours = (cc(100)))

p1+p2



sum(tcga_obj_sub$NEU_xcell > median(tcga_obj_sub$NEU_xcell) &
      tcga_obj_sub$NEU_mcp > median(tcga_obj_sub$NEU_mcp) &
      tcga_obj_sub$NEU_quantiseq > median(tcga_obj_sub$NEU_quantiseq)
)
sum(tcga_obj_sub$NEU_xcell < median(tcga_obj_sub$NEU_xcell) &
      tcga_obj_sub$NEU_mcp < median(tcga_obj_sub$NEU_mcp) &
      tcga_obj_sub$NEU_quantiseq < median(tcga_obj_sub$NEU_quantiseq)
)

tcga_obj_sub$NEU = "1hetero"
tcga_obj_sub$NEU[tcga_obj_sub$NEU_xcell > quantile(tcga_obj_sub$NEU_xcell, 0.3333) & tcga_obj_sub$NEU_mcp > quantile(tcga_obj_sub$NEU_mcp, 0.3333) & tcga_obj_sub$NEU_quantiseq > quantile(tcga_obj_sub$NEU_quantiseq, 0.3333)] = "2High"
tcga_obj_sub$NEU[tcga_obj_sub$NEU_xcell < quantile(tcga_obj_sub$NEU_xcell, 0.6666) & tcga_obj_sub$NEU_mcp < quantile(tcga_obj_sub$NEU_mcp, 0.6666) & tcga_obj_sub$NEU_quantiseq < quantile(tcga_obj_sub$NEU_quantiseq, 0.6666)] = "0Low"
table(tcga_obj_sub$NEU)

p3 = DimPlot(tcga_obj_sub, reduction = "tsne", label = T, pt.size = .5, 
             group.by = "NEU",  raster = F)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
tcga_obj_sub$NEU_xcell_norm = range01(tcga_obj_sub$NEU_xcell)
tcga_obj_sub$NEU_mcp_norm = range01(tcga_obj_sub$NEU_mcp)
tcga_obj_sub$NEU_quantiseq_norm = range01(tcga_obj_sub$NEU_quantiseq)

quantile(tcga_obj_sub$NEU_xcell_norm)
quantile(tcga_obj_sub$NEU_mcp_norm)
quantile(tcga_obj_sub$NEU_quantiseq_norm)


tcga_obj_sub$NEU_consensus = rowMeans(cbind(tcga_obj_sub$NEU_xcell_norm, tcga_obj_sub$NEU_mcp_norm, tcga_obj_sub$NEU_quantiseq_norm))
quantile(tcga_obj_sub$NEU_consensus)


p4 = FeaturePlot(tcga_obj_sub, features = c("NEU_consensus"), pt.size = .01, ncol = 1, reduction = "tsne") & theme_bw() & theme(aspect.ratio=1)

cc <- colorRampPalette(c("#352a86", "#095cd8", "#46b896", "#e7ba4a", "#f8fa0d"))
p5 = Nebulosa::plot_density(tcga_obj_sub, "NEU_consensus", size = 0.5, reduction = "tsne") & theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank()) &
  scale_color_gradientn(colours = (cc(100)))

p1+p2+p3+p4+p5









############
# Single-cell visualization
DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .000000000000000001, 
        group.by = "celltype_l1", label.size = 0, raster = F)

All_obj_time1$type2.1 = "Other"; All_obj_time1$type2.1[All_obj_time1$type2 == "Healthy/Adjacent" & All_obj_time1$site != "PBMC"]="Adjacent/Healthy"
table(All_obj_time1$type2.1)
All_obj_time1$type2.2 = "Other"; All_obj_time1$type2.2[All_obj_time1$type2 == "Cancer"]="Cancer"
All_obj_time1$type2.3 = "Other"; All_obj_time1$type2.3[All_obj_time1$type2 == "Metastasis"]="Metastasis"
All_obj_time1$type2.4 = "Other"; All_obj_time1$type2.4[All_obj_time1$site == "Circulating"]="Circulating"

DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .000000000000000001, group.by = "type2.2", 
        cols = c("grey30", "grey97"), label.size = 0, raster = F) + 
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .000000000000000001, group.by = "type2.1", 
          cols = c("grey45", "grey97"), label.size = 0, raster = F) + 
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .000000000000000001, group.by = "type2.3", 
          cols = c("grey45", "grey97"), label.size = 0, raster = F) + 
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  DimPlot(All_obj_time1, reduction = "umap", label = TRUE, pt.size = .000000000000000001, group.by = "type2.4", 
          cols = c("grey60", "grey97"), label.size = 0, raster = F) + 
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())





#barplot
dittoBarPlot(All_obj_time1, "celltype_l1", group.by = "type", scale = "percent")





#############
# published neutrophil states


library(UCell)
signature = list(
  Ecotyper_S01 = c("AZU1","PCSK5","KIR3DX1","CCDC91","PTP4A2","SCGB1A1","ZNF121","DGAT2","NRM","GTF2F2"),
  Ecotyper_S02 = c("CCL23","HBA2","QKI","HBB","SFTPC","HBA1","EDNRB","PTX3","CDS1","S1PR1"),
  Ecotyper_S03 = c("SIGLEC5","CLEC12A","USP10","SIRPB1","FPR1","APMAP","RGL4","SH3KBP1","MNDA","CR1"),
  hNeutro1 = c("PADI4","MMP9","S100A12","HBB","ARG1","FCN1","CKAP4","CDA","HMGB2","VIM","PLBD1","APMAP","S100P","RGL4","IL18RAP","AC245128.3","AGFG1","ADD3","CLEC4D","CRISP3","CALM1","MARC1","RFLNB","LINC02207","C1orf162","S100A8","PGLYRP1","KIF1B","TFDP1","FOSL1","VNN1","TSPO","ID3","ADAM19","PRKAR2A","PGD","PSTPIP2","FBXW2","OSBPL9","AC068473.5","S100A4","HIPK2","PRNP","GYG1","AL034397.3","SLC36A4","PADI2","EIF4G2","HK3","RBP7","TMCC3","CRISPLD2","NR4A3","CLC","DBNL","RBL2","MANSC1","BMX","TUBB4B","SYTL3","PFKFB2","MAEA","PPP2R5A","S100A9","DHCR7","STK38L","CLEC12A","YWHAB","SDF2","CTDSP1","RAB27A","SBNO2","S100A6","MAST3","WASF2","AL021707.6","MAP3K5","RSBN1","PELO","LAMTOR5","SORT1","CAP1","ATM","IL1R2","STK38","SYNE2","FKBP9","LINC00694","SEPT9","ATP11B","STK17B","TLE3","C19orf38","TUBA4A","CTDP1","MYL6","DOK3","ROCK1","PINK1","ACSL1"),
  hNeutro2 = c("ISG15","RSAD2","IFIT2","IFIT3","IFIT1","XAF1","HERC5","OASL","MX1","PARP14","IFI6","UBE2L6","PML","RNF213","IRF7","IFI44L","GBP5","OAS2","GBP4","TMEM123","ANKRD22","STAT2","EIF2AK2","SAMD9","TNFSF10","STAT1","OAS3","SAMD9L","DDX60L","SERPING1","TRIM22","IFIT5","IFIH1","EPSTI1","GBP1","NT5C3A","BATF2","NUB1","IFITM3","CMPK2","BST2","ERAP2","SAMHD1","PARP9","APOL6","DDX58","TAGAP","SPTA1","TRIM21","PHF11","IRF1","IFI44","GLRX","PARP12","C5orf56","ZBP1","AC087672.3","HES4","LYSMD2","U62317.5","LACTB","LY6E","TNFSF13B","RBBP6","PSMB9","ZNFX1","HSH2D","NCOA7","MICB","TRAFD1","IFI16","ATF5","IFI35","HERC6","DNAJB12","FFAR2","C8orf49","ARFGEF2","TAP1","SCARB2","IFITM1","REC8","DHX58","CNOT8","TRAPPC13","LINC01840","NEURL3","PRELID2","ACOT9","DDX60","SP110","NPFFR1","MRPL30","ABCB9","FEZ2","PLSCR1","MT2A","XRN1","GLB1L","GEMIN8"),
  hNeutro3 = c("ARHGEF40","CDC42EP2","CPPED1","CXCR2","MPZL1","CFD","SNUPN","CASC3","AMPD2","VMP1","NCOA4","JAML","RASSF3","MICAL2","FRAT2","NCOA1","EGLN1","CCR4","GLT1D1","FAM129A","MPPE1","PIGX","OAZ2","BTN3A1","SELPLG","LINC01506","ALPL","CASS4","IDS","NDE1","RN7SL455P","VPS35","SH2B2","RALBP1","PECAM1","HOTAIRM1","LINC01460","TAGLN2","STIM1","BACE1","PPM1F","HIST1H2BC","ALAD","KIAA2013","AKNA","TSC22D3","LRRC4","TPR","AC012368.1","LARP1B","LINC02555","AC090921.1","ITM2B","TUBA1B","ITPK1","TNRC6B","EIF1AY","SNX6","KLRG1","CENPF","ARNTL","POLB","R3HDM4","DPYSL5","NUDT16","AC012645.4","ZC3H7B","ST3GAL1","ZNF236","CAPN15","PTPRC","NSFL1C","RUNX2","LINC02242","SDE2","FAM212B","HLA-DQB1","JMJD1C","SPATA13","KRT23","GNG10","FOXN3","FRY","GKAP1","TLR1","AC026401.3","DHRS12","COL18A1","CYBRD1","PPIA","FAM157C","NCAM1","THG1L","MYO1F","TXNIP","ST8SIA4","ST3GAL4","SLC8A1","CD302","CYTH4"),
  hNeutro4 = c("HLA-DRA","DDX3Y","ARRB1","HLA-DRB1","NR4A2","RHOH","CTSC","ATP1B3","RBM3","SERPINB9","ARL6IP5","CALHM6","PPIF","SIGLEC10","RNF141","ISG20L2","ATP1A1","KMT5A","CAND1","SPIDR","CLEC4A","IVNS1ABP","CCL4L2","ETV5","GK-IT1","CD74","CYTOR","IL1RAP","ZBTB8OS","AC016831.1","ROBO1","SLC7A5","LINC02078","FCGR2B","LCP2","CLEC7A","TRAF3IP2","TUFM","PYCR2","AC068473.3","SCAMP1-AS1","EIF4A3","AL136295.7","TNFAIP3","UBE2W","BIRC3","HIC1","SUN2","MIR4435-2HG","LTBR","IRF2BP2","AC020916.1","GREM1","SPRTN","ITGAX","MEF2C","EPHA1","NGLY1","CDC42EP3","CMTM6","DNAJC24","RPS15","FAM126A","JUNB","RNASEH2B","CYBB","EHD1","PPT1","MPEG1","CEBPB","ARHGAP18","TMEM209","SPAG9","NFKBIA","ADGRE2","CD69","GZF1","BRDT","SBDS","SPP1","HDGF","Metazoa_SRP","SCARF1","IL1B","LST1","NFKBIE","USP36"),
  hNeutro5 = c("PI3","SPINK1","MT1G","FNIP2","IL3RA","CSTB","LGALS3","CAST","PLPP3","SLPI","AC015912.3","TPRA1","GPR84","CSF1","CCL3","SNAPC1","TNFSF15","RMND5A","CCDC93","GNPDA1","ATP6V1C1","NPC1","RHEB","GBE1","TGM3","APLP2","TGM2","DARS","ATP6V1F","USF2","TMEM251","CNST","UQCRC2","BNIP3L","SEC22B","TRGC2","CCL20","IRAK1","ZNF438","ID2","JUN","C12orf49","TPI1","TRIM33","TBC1D7","CCL2","CDCP1","IRAK2","CORO1C","AL031316.1","RPS28","PKM","GLMP","ATP6V1G1","MIF","RRAGD","SGSH","SSR2","GSTO1","PLEKHB2","TRPM2","DENND5A","GABARAPL2","TXN","OTOA","TMED10","MGAT1","MFSD2A","MGLL","HMOX1","GAPDH","VPS18","BRI3","ATP6AP1","PDXK","SLC43A3","ZNF316","FKBP15","NPL","MIR222HG","KCNAB2","STRBP","ANKRD12","RALGDS","KRCC1","FNDC3A","ZFAND5","HSPA1A","ACVR1B","OSTM1","SLC11A2","MAFG","KLHDC8B","YBX1","CD48","DNAJB1","MAPK6","PHLDA1","AMDHD2","CD83"),
  CC_TAN1 = c("IL1RN", "CD44", "RHOH", "ELOC", "RIPK2"),
  CC_TAN2 = c("HLA-DRA", "HLA-DRB1", "CD74", "SNRPG", "HLA-DMB"),
  CC_TAN3 = c("ASAH1", "CSTB", "LGALS3", "GNPDA1", "PLEKHB2"),
  CC_TAN4 = c("ATP5MC2", "RPL23", "RPN2", "RPL3", "RPS12")
)
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(All_obj_time1@assays$RNA@data, features = signature, name = "")))

#signature.matrix.neu is the signature score matrix



