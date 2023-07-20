suppressPackageStartupMessages(library("reticulate"))
use_condaenv("py38")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggplot2))

df_full=data.table::fread("./df_full_compare_new.csv.gz")
df_full$paper=plyr::mapvalues(df_full$Group,c("our","nature","cir","Science","Cui"),to=c("This study","Litviňuková et al","Tucker et al","Reichart et al","Cui et al"))
df_full=as.data.frame(df_full)
df_full$paper=factor(df_full$paper,levels = c("Tucker et al","Cui et al","Litviňuková et al","Reichart et al","This study"))
table(df_full$paper)

get_plot=function(y_var="n_genes_by_counts",ylab_title="nGene",ggtitle0=NULL){
  df_full$y=df_full[,y_var]
  p1=ggplot(data=df_full,aes(x=paper,y=y,fill=paper))+
    geom_violin(scale = "width")+
      scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Set2"))+
    geom_boxplot(width=0.15,color="blue",outlier.color = NA,size=0.2)+
    geom_point(data=df_full%>%group_by(paper)%>%summarise(mean=mean(y)),aes(x=paper,y=mean),size=2,color="white")+
  #p1=ggviolin(df_full, "paper", "y",  color = "paper",
   #palette = RColorBrewer::brewer.pal(5, "Set2"),
   #add = "boxplot")+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 30,hjust=1,size=15),
          axis.text.y=element_text(size=15),
          axis.title.y = element_text(size=16),
          strip.text = element_text(size=20,color="black"),
          )+xlab("")+
    ylab(ylab_title)+
    ggtitle(ifelse(is.null(ggtitle0),ylab_title,ggtitle0))+
    theme(plot.title = element_text(size=20,hjust=0.5))
  
    return(p1)
}
comparisons=list(c("This study","Tucker et al"),
                 c("This study","Cui et al"),
                 c("This study","Litviňuková et al"),
                 c("This study","Reichart et al"))
p1=get_plot("n_genes_by_counts",'nGenes')+ggpubr::stat_compare_means(comparisons =comparisons)
p2=get_plot(y_var="total_counts",'nCounts')+ggpubr::stat_compare_means(comparisons =comparisons)
p3=get_plot(y_var="pct_counts_tf",'tf.percent(%)',"tf.percent")+scale_y_continuous(trans = "sqrt")+ggpubr::stat_compare_means(comparisons =comparisons)
p4=get_plot(y_var="pct_counts_mt",'mt.percent(%)','mt.percent')+scale_y_continuous(trans = "sqrt")+ggpubr::stat_compare_means(comparisons =comparisons)
p5=get_plot(y_var="pct_counts_rib",'rib.percent(%)','rib.percent')+scale_y_continuous(trans = "sqrt")+ggpubr::stat_compare_means(comparisons =comparisons)
p6=get_plot(y_var="total_counts_tf",'nTF_Counts','nTF_Counts')+scale_y_continuous(trans = "sqrt")+ggpubr::stat_compare_means(comparisons =comparisons)
p=patchwork::wrap_plots(list(p1,p2,p3,p4,p5,p6),ncol=3)
ggsave("./paper_output/revised_figures/comparison.pdf",p,width=18,height=10)
df1=df_full%>%group_by(paper)%>%summarise(mean_nGene=mean(n_genes_by_counts),median_nGene=median(n_genes_by_counts))
df2=df_full%>%group_by(paper)%>%summarise(mean_nCount=mean(total_counts),median_nCount=median(total_counts))
df3=df_full%>%group_by(paper)%>%summarise(mean_tf=mean(pct_counts_tf),median_tf=median(pct_counts_tf))
df4=df_full%>%group_by(paper)%>%summarise(mean_mt=mean(pct_counts_mt),median_mt=median(pct_counts_mt))
df5=df_full%>%group_by(paper)%>%summarise(mean_rib=mean(pct_counts_rib),median_rib=median(pct_counts_rib))
df6=df_full%>%group_by(paper)%>%summarise(mean_tf_count=mean(total_counts_tf),median_tf_count=median(total_counts_tf))
df0=plyr::join_all(list(df1,df2,df3,df4,df5,df6), by='paper', type='left')
openxlsx::write.xlsx(df0,"./paper_output/revised_figures/comparison_related.xlsx")
