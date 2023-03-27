
##########
#bar plot visualization


propmat_sub = read.table("input.txt", header=T, sep="\t")

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

propmat_sub2 = data_summary(propmat_sub, varname="Statistic", 
                            groupnames=c("Group"))


library(ggpubr)
library(ggplot2)
ggplot(data = propmat_sub, aes(x = Group, y = Statistic))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(data = propmat_sub, aes(fill = Group, color = Group), size=3, #, shape = group2
             position = position_jitterdodge(dodge.width = 0.2), 
             alpha=0.7)+
  theme_bw()+
  theme(aspect.ratio=1)+
  theme(legend.position="right", legend.box = "right", axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_compare_means(comparisons = list(
    c("group1","group2")
  ),method = "t.test",label = "p.signif") +
  NULL



ggplot(data = propmat_sub, aes(x = Group, y = Statistic))+
  geom_bar(data = propmat_sub2, aes(x = Group, y = Statistic, fill = Group),
           stat="identity", color="black", fill = "grey95", width = .45,
           position=position_dodge()) +
  geom_errorbar(data = propmat_sub2, aes(ymin=Statistic, ymax=Statistic+sd), width=.4,
                position=position_dodge(.9)) +
  
  geom_point(data = propmat_sub, aes(fill = Group, color = Group), size=3, color = "black", shape = 21,
             position = position_jitterdodge(dodge.width = 0.2), 
             alpha=1)+
  coord_cartesian(ylim = c(1000, 5000)) +
  theme_bw()+
  theme(aspect.ratio=2)+
  theme(legend.position="none", legend.box = "none", axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c(c(my36colors))) +
  stat_compare_means(comparisons = list( 
    c("group1","group2")
  ),method = "t.test",label = "p.signif") + #,label = "p.label"
  NULL

