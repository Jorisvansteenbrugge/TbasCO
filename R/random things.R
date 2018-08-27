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


less8 <- which((apply(Model_Module_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)<=7) == T) %>% names
Model_Module_List_8 <- Model_Module(RNAseq.data, trait.attributes, 39, less8, Bin_Order_Index, Yrange, bkgd.traits)

a <- which(Model_Module_List$Model_Sig_Matrix['16',] == T) %>% names
b <- which((apply(Model_Module_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)<=3) == T) %>% names

Module_Names_nd<-intersect(a,b)

A <- which(Model_Module_List$Model_Sig_Matrix['16',] == T) %>% names
B <- which(apply(Model_Module_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)>3 &
              (apply(Model_Module_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)<=7) == T) %>% names

Module_Names_mid<-intersect(A,B)

Aa <- which(Module_Model_List$Model_Sig_Matrix['16',] == T) %>% names
Bb <- which((apply(Module_Model_List$Model_Sig_Matrix,2,sum,na.rm=TRUE)>=3) == T) %>% names

Module_Names_eco<-intersect(Aa,Bb)


Module_Names_Polymer_Metabolism <- c("M00001","M00307", "M00579","PHA", "M00172")
Module_Names_Purine_Metabolism <-c("M00048","M00049","M00050","M00546")
Module_Names_Pyrimidine_Metabolism <-c("M00051","M00052","M00053","M00046")
Module_Names_PP_Metabolism <- c(Module_Names_Purine_Metabolism,Module_Names_Pyrimidine_Metabolism)
Module_Names_ABC_transporters<-c("M00185","M00186","M00189","M00190","M00191","M00194","M00196","M00198","M00199","M00200","M00201","M00202","M00203","M00204","M00205","M00206","M00208","M00209","M00210","M00212","M00213","M00214","M00215","M00216","M00217","M00218","M00219","M00220","M00222","M00223","M00224","M00225","M00226","M00227","M00228","M00229","M00230","M00231","M00232","M00233","M00234","M00235","M00237","M00238","M00239","M00240","M00241","M00242","M00245","M00246","M00249","M00250","M00251","M00252","M00253","M00255","M00256","M00257","M00259","M00298","M00299","M00300","M00301","M00302","M00314","M00316","M00317","M00318","M00319","M00320","M00321","M00322","M00323","M00324","M00325","M00348","M00349","M00423","M00435","M00436","M00437","M00438","M00439","M00440","M00442","M00491","M00566","M00581","M00582","M00583","M00584","M00585","M00586","M00587","M00589","M00590","M00591","M00592","M00593","M00599","M00600","M00601","M00602","M00603","M00604","M00605","M00606","M00607","M00619","M00634","M00635","M00731","M00732","M00739","M00747","M00762","M00791","M00792","M00813","M00817")
Module_Names_Nitrogen_Metabolism <- c("M00175","M00531","M00530","M00529","M00528","M00804")
Module_Names_CC_Metabolism <- c("M00001","M00002","M00003","M00307","M00009","M00010","M00011","M00004","M00006","M00007","M00580","M00005","M00008","M00308","M00633","M00309")
Module_Names_ATP_synthesis <- c("M00144","M00143","M00146","M00147","M00149","M00150","M00148","M00162","M00151","M00152","M00154","M00155","M00153","M00417","M00416","M00156","M00157","M00159")
Module_Names_Phosphotransferase_system<-c("M00265","M00809","M00267","M00266","M00806","M00282","M00269","M00271","M00272","M00270","M00303","M00268","M00273","M00306","M00274","M00305","M00281","M00275","M00280","M00279","M00807","M00276","M00764","M00304","M00278","M00277","M00287","M00610","M00283")
Module_Names_MOI_transport_system <- c("M00185","M00189","M00423","M00186","M00438","M00321","M00188","M00435","M00436","M00437","M00190","M00191","M00299","M00300","M00193","M00301","M00302","M00208","M00209","M00442","M00192")


Model_Module_List <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names, Bin_Order_Index, Yrange, bkgd.traits)
Model_Module_List_nd <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_nd, Bin_Order_Index, Yrange, bkgd.traits)
Model_Module_List_mid <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names, Bin_Order_Index, Yrange, bkgd.traits)
Model_Module_List_high <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names, Bin_Order_Index, Yrange, bkgd.traits)

Model_Module_List_All <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_All, 39, Module_Names)

Model_Module_List_PP_Metabolism <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_PP_Metabolism, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_PP_Metabolism, 39, Module_Names_PP_Metabolism)


Model_Module_List_All <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_All, 39, Module_Names)


Model_Module_List_Phosphotransferase_system <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_Phosphotransferase_system, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_Phosphotransferase_system, 39, Module_Names_Phosphotransferase_system)

Model_Module_List_MOI_transport_system <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_MOI_transport_system, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_MOI_transport_system, 39, Module_Names_MOI_transport_system)


Model_Module_List_ABC_transporters <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_ABC_transporters, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_ABC_transporters, 39, Module_Names_ABC_transporters)

Model_Module_List_CC_Metabolism <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_CC_Metabolism, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_CC_Metabolism, 39, Module_Names_CC_Metabolism)


Model_Module_List_ATP_synthesis <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_ATP_synthesis, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_ATP_synthesis, 39, Module_Names_ATP_synthesis)

Model_Module_List_Polymer_Metabolism <- Model_Module(RNAseq.data, trait.attributes, 39, Module_Names_Polymer_Metabolism, Bin_Order_Index, Yrange, bkgd.traits)
Plot_Model_Module(Model_Module_List_Polymer_Metabolism, 39, Module_Names_Polymer_Metabolism)


sapply(sapply(sub_modules[], length), max)

