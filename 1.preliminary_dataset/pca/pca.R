setwd("~/Desktop/PanPan/POPGEN_short-read_Illumina/PCA/")
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(RcppCNPy)
plot_pca1<-function(matrix,components, n1, n2, e1,e2){
  names<-read.delim("allchimps_bams.txt", header = FALSE)
  dataset<-read.delim2('Coriell_dataset.txt', sep="\t", header=TRUE)
  names$CODE <- gsub("/global/scratch/users/joana_rocha/PANPAN/CHIMPS_lowcov/bam_final/", "", names$V1)
  names$CODE <- gsub(".sorted.bam", "", names$CODE)
  covmatrix<-read.delim(matrix, sep = " ",header=FALSE)
  e<-eigen(covmatrix)
  print(e$values[e1]/sum(e$values)*100)
  print(e$values[e2]/sum(e$values)*100)
  pca <- as.data.frame(e$vectors[,components])
  vars<-apply(pca, 2, var)
  vars/sum(vars)
  colnames(pca) <- c("x", "y")
  pca$CODE<-names$CODE
  ingroup_pca <- merge(dataset,pca,by="CODE")
  ingroup_pca$mtDNASubspecies<-factor(ingroup_pca$mtDNASubspecies, levels = c("P. troglodytes troglodytes", "P. troglodytes verus", "P. troglodytes schweinfurthii"))
  ggplot(ingroup_pca, aes(x=x, y=y, colour=mtDNASubspecies), main=NULL) + geom_point(size=0.8) + 
  xlab(n1) + ylab(n2) + theme_test(base_size = 11) + scale_color_manual(values = colors1) 
  }
colors1 <- c('#4c5d4c','#9dced9', '#ffb35a')
# '#4c5d4c' -> central african (Pantrog trog)
# '#9dced9' -> western African (Pantrog verus)
# '#fe604c' -> nigeria cameroon (Pantrog elioti)
# '#ffb35a' -> eastern (Pantrog schweinfurthii)

myplot1<-plot_pca1('allchimps_snps_fromBam.cov',  1:2,"PC1", "PC2", 1,2) +
  geom_label_repel( aes(label = CODE),box.padding=0.15, point.padding = 0.01, label.size = 0.001, segment.size =0.5, size = 5) 
myplot2<-plot_pca1('allchimps_snps_fromBam_minmaf0.05.cov',  1:2,"PC1", "PC2", 1,2) +
  geom_label_repel( aes(label = CODE),box.padding=0.15, point.padding = 0.01, label.size = 0.001, segment.size =0.2, size = 3) 

#my_components<-ggarrange(myplot1, myplot2, nrow=1, common.legend = TRUE, legend="bottom", labels = c("A"),  font.label = list(size = 16))
#my_components_vr<-ggarrange(myplot3, myplot4, nrow=1, common.legend = TRUE, legend="bottom", labels = c("B"), font.label = list(size = 16))
#print(my_components_vv)
#print(my_components_vr)
#final<-ggarrange(my_components_vv,my_components_vr, nrow = 2)
#+ theme(legend.title = element_text(color="black", size=10), legend.text=element_text(size = 10))

ggsave("myplot.tiff", myplot1, width = 5, height =3, dpi=800)

#mainfig<-ggarrange(myplot1, legend="bottom", labels = c("A", "B"), ncol= 2,font.label = list(size = 12))
#ggsave("mainfig.tiff", mainfig, "tiff", dpi=300, width=11, height=6, units="in")

