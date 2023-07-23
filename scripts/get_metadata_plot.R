
celltype_color_defined=c("CM"="#c10023",
                         "FB"="#FE0092",
                          "FB1"="#FE0092",
                         "FB2"="#FF66CF",
                         "EC"="#008e17",
                         "EC1"="#008e17",
                         "EC2"="#627B39",
                         "EC3"="#B6CF6A",
                         "endoEC"="#00cc99",
                         "SMC"="#A55296",
                         "Peri"="#E8979C",
                         "Peri1"="#E8979C",
                         "Peri2"="#FF7947",
                         "Peri3"="#EFB8FF",
                         "Myeloid"="#0099cc",
                         "Lymphoid"="#0053c8",
                         "Neu"="#B44B4D",
                         "Doublet"="grey")#"Neutrophils"="#bc9000"

obj=readRDS("sampleid_specify.rds")#
Idents(obj)="maincelltype"
celltype_colors_used=celltype_color_defined[levels(Idents(obj))]

get_plot=function(obj,color.used,group.by="maincelltype",label_cellnumber=T){	
  colors_tmp=color.used
  Idents(obj)=group.by
  df0=cbind(as.data.frame(obj[["tsne"]]@cell.embeddings),obj@meta.data)
  df0$group0=df0[,group.by]
  df0_median=df0%>%group_by(group0)%>%summarise(x=median(tSNE_1),y=median(tSNE_2),n=n())
  df0_median$color=plyr::mapvalues(df0_median$group0,names(colors_tmp),colors_tmp)
  if(label_cellnumber){
    df0_median$label= paste0(df0_median$group0," (",scales::comma(as.numeric(df0_median$n),accuracy = 1),')')
  }else{
    df0_median$label=df0_median$group0
  }
  p_celltype=ggplot(df0) + 
    geom_point(aes_string(x="tSNE_1", y="tSNE_2", color="group0"),size=0.1) +
  scale_color_manual(values=colors_tmp)+
  guides(color=guide_legend(override.aes = list(size=5),ncol=1))+
  theme(legend.text = element_text(hjust=0))+
  geom_point(data=df0_median,aes(x=x,y=y),size=3,color="lightgrey",alpha=0.3)+
  guides(color=guide_legend(override.aes = list(size=5),ncol=1))+
  ggrepel::geom_label_repel(data=df0_median,aes(x=x,y=y,label=label),
                            #fill=alpha(df0_median$color,0.3),
                            fill=df0_median$color,
                            color="black",
                            segment.colour = "black",
                            fontface = 'bold',
                            point.padding = unit(0.5, "lines"),
                   label.size = NA, 
                   size=5,
                   alpha = 0.15, 
                   label.padding = 0.2,
                   na.rm=TRUE,seed=12)+ #geom_text(data=df0_median,aes(x=x,y=y,label=celltype),size=5)+
    ggrepel::geom_label_repel(data=df0_median,aes(x=x,y=y,label=label),
                            #fill=alpha(df0_median$color,0.3),
                            fill=NA,
                            color="black",
                            segment.colour = "black",
                            fontface = 'bold',
                            point.padding = unit(0.5, "lines"),
                   label.size = NA, 
                   size=5,
                   alpha = 0.8, 
                   label.padding = 0.2,
                   na.rm=TRUE,seed=12)+
    cowplot::theme_cowplot()+
  theme(legend.text=element_text(hjust=0),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.01, 'cm'))+
    theme_00+
    xlab("")+
  ylab("")
  return(p_celltype)
}
p=get_plot(obj,color.used = celltype_colors_used,group.by = "maincelltype")
p
