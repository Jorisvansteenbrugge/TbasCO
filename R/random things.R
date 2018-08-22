library(gplots)
heatmap.2(sbs.trait.attributes[,as.numeric(which(colSums(sbs.trait.attributes[,])>3))])
heatmap.2(pa_trait[which(rowSums(pa_trait) > 3),])



hclust(dist(sbs.trait.attributes[,as.numeric(which(colSums(sbs.trait.attributes[,])>3))])) %>% plot(main = "Attributes")

# Traits
pa_trait <- apply(RNAseq.data$features$trait_presence_absence, 2, as.numeric)
pa.clust <- pa_trait[which(rowSums(pa_trait) > 3),]%>% t %>% dist %>% hclust
plot(pa.clust, labels = genome.taxonomy.phylum)

pa.clust.attributes <- sbs.trait.attributes[,as.numeric(which(colSums(sbs.trait.attributes[,])>3))] %>% dist %>% hclust

plot(pa.clust.attributes, labels = genome.taxonomy.phylum)



# overlap in niche specific 16 & 39
# ta is from the Plot_Redundancy_Traits function in utility.R
which(RNAseq.data$features$trait_presence_absence[names(which(ta.pa<4)), "16"] ==1)[which(RNAseq.data$features$trait_presence_absence[names(which(ta.pa<4)), "16"] ==1)%in%which(RNAseq.data$features$trait_presence_absence[names(which(ta.pa<4)), "39"] ==1)]




Module_Names <- RNAseq.data$features$trait_presence_absence[,'39'] %>% which(. == T) %>% names




Figure_X<niche - barplot(sort(apply(Module_Model_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)), xaxt ='n')
labels <- colnames(Module_Model_List$Model_Comparison_Matrix)[Module_Model_List$Module_Order_Index]
text (cex=0.75, x=Figure_X-.25, y=-.6, labels, xpd=TRUE, srt=90)


a <- which(Module_Model_List$Model_Sig_Matrix['16',] == T) %>% names
b <- which((apply(Module_Model_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)<=3) == T) %>% names

intersect(a,b)

A <- which(Module_Model_List$Model_Sig_Matrix['16',] == T) %>% names
B <- which(apply(Module_Model_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)>3&
              (apply(Module_Model_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)<=10) == T) %>% names

intersect(A,B)
