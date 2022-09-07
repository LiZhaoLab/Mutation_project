setwd("~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/")
library(Seurat)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(EnhancedVolcano)
library(cowplot)
library(corrplot)


Old1<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/Other data/Old1_sc/filtered_feature_bc_matrix")
Old1 <- CreateSeuratObject(Old1,min.cells=3, min.features =200)
Old1$rep<-"Old1"


Old3a<-Read10X("Old3a_run/")
Old3a<-CreateSeuratObject(Old3a,min.cells=3, min.features =200)
Old3a$rep<-"Old3a"


Old3b<-Read10X("Old3b_run/")
Old3b<-CreateSeuratObject(Old3b,min.cells=3, min.features =200)
Old3b$rep<-"Old3b"



Rep1<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/R517_3reps/Evanlib1_run/")
Rep1 <- CreateSeuratObject(Rep1,min.cells = 3, min.features =200)
Rep1$rep<-"Rep1"


Rep2<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/R517_3reps/Evanlib2_run/")
Rep2 <- CreateSeuratObject(Rep2, min.cells = 3,min.features =200)
Rep2$rep<-"Rep2"


Younglib3<-Read10X("YoungLib3_run/")
Younglib3<-CreateSeuratObject(Younglib3,min.cells=3, min.features =200)
Younglib3$rep<-"Younglib3"



#normalize and identify variablefeatures
testis.list<-lapply(X = list(Old1, Old3a,Old3b, Rep1, Rep2, Younglib3), FUN = function(x){
 # x<-SCTransform(x)
  #x<-RunPCA(x)
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, nfeatures = 2000)
  x<-ScaleData(x)
  x<-RunPCA(x)
 # x<-RunUMAP(x, dims = 1:10)
  x<-SCTransform(x)
})



features <- SelectIntegrationFeatures(object.list = testis.list)
anchors <- FindIntegrationAnchors(object.list = testis.list, anchor.features = features)
# this command creates an 'integrated' data assay
testis.list <- IntegrateData(anchorset = anchors)

testis.list<-SCTransform(testis.list )
testis.list <- FindVariableFeatures(testis.list, assay = "SCT")

testis.list <- RunPCA(testis.list,verbose = FALSE, npcs = 15, assay = "SCT")
#testis.list<-RunCCA(testis.list)
testis.list <- RunUMAP(testis.list, reduction = "pca", assay = "SCT", dims = 1:11)

testis.list<-FindNeighbors(testis.list, reduction = "pca", assay = "SCT" )
testis.list<-FindClusters(testis.list, resolution = 1.2 )
#FeaturePlot(testis.list, split.by = "rep",features   = c("nCount_RNA","nFeature_RNA"))
DimPlot(testis.list, split.by = "Age", label = T)






#saveRDS(testis.list, "~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps.RDS")
FeaturePlot(testis.list, split.by = "rep",features = c("p-cup","aub","soti","fzo"))

FeaturePlot(testis.list, split.by = "rep",features = c("p-cup"))
FeaturePlot(testis.list, split.by = "rep",features = c("aub"))
FeaturePlot(testis.list, split.by = "rep",features = c("bam"))
FeaturePlot(testis.list, split.by = "rep",features = c("dlg1"))

FeaturePlot(testis.list, split.by = "rep",features = c("soti"))
FeaturePlot(testis.list, split.by = "rep",features = c("Dpy-30L2"))
FeaturePlot(testis.list, split.by = "rep",features = c("Prosalpha6"))
FeaturePlot(testis.list, split.by = "Age",features = c("MtnA"))

Idents(testis.list)<-factor(Idents(testis.list), levels = c("Hub cells","Cyst cells","Epithelial cells","Accessory gland","GSC, Early spermatogonia","Late spermatogonia","Early spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"))
testis.list$Age<-ifelse(testis.list$rep %in% c("Old1","Old3a","Old3b"),"Old","Young")

DimPlot(testis.list, split.by = "Age", label = F)
DotPlot(testis.list, features = c("p-cup","aub","Prosalpha6T", "His2Av", "ProtB", "Dpy-30L2", "fzo","twe","dlg1","MtnA", "Sems", "Fas3", "Rab11"))
DotPlot(testis.list,split.by = "Age", features = c("p-cup","soti","fzo","twe"))
DotPlot(testis.list, features = c("p-cup","soti","fzo","twe"))



new.idents<-c("Epithelial cells",#0
              "Early spermatids",#1
              "Early spermatocytes",#2
              "Late spermatogonia",#3
              "Epithelial cells", #4
              "Late spermatocytes",#5
              "Epithelial cells",#6
              "Epithelial cells",#7
              "Early spermatids",#8
              "Cyst cells",#9
              "Late spermatocytes",#10
              "Early spermatids",#11
              "Hub cells",#12
              "Epithelial cells",#13
              "Early spermatids",#14
              "Cyst cells",#15
              "Epithelial cells",#16
              "Late spermatogonia",#17
              "Late spermatocytes",#18
              "Hub cells",#19
              "Late spermatogonia",#20
              "Late spermatocytes",#21
              "Epithelial cells",#22
              "Late spermatogonia",#23
              "GSC, Early spermatogonia",#24
              "Cyst cells",#25
              "Early spermatids",#26
              "Early spermatids",#27
              "Epithelial cells",#28
              "Accessory gland",#29
              "Late spermatids",#30
              "Early spermatids",#31
              "Late spermatogonia"#32
              )



names(new.idents) <- levels(testis.list)
testis.list<-RenameIdents(object= testis.list, new.idents )
testis.list$ident<-Idents(testis.list)

DimPlot(testis.list, split.by = "Age")
DotPlot(testis.list, features = c("p-cup","aub","Prosalpha6T", "His2Av", "ProtB", "Dpy-30L2", "fzo","twe","dlg1","MtnA", "Sems", "Fas3", "Rab11"))




#saveRDS(testis.list, "~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps.RDS")
testis.list<-readRDS("~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps.RDS")
testis.list$ident<-Idents(testis.list)


#maybe try a loop to get tidier data
Idents(testis.list)<-"ident"
par(mfrow=c(2,3), xpd = T)
#pdf("Supplemental_Figure_1.pdf", width=8, height = 8)
for (i in c("GSC, Early spermatogonia", "Late spermatogonia","Early spermatocytes","Late spermatocytes","Early spermatids","Late spermatids")){
  Avgexpsplit <-AverageExpression(subset(testis.list, ident ==i), group.by = c("rep"))$RNA
  cormat<-cor(Avgexpsplit)
  colnames(cormat)<-c("Old1","Old2","Old3","Young1","Young2","Young3")
  rownames(cormat)<-c("Old1","Old2","Old3","Young1","Young2","Young3")
  corrplot(cormat, method = 'number',tl.col = "black", title = i, mar=c(0,0,1,0) )
  
}
#dev.off()



#Figure 1

DimPlot(testis.list, label = T)
Idents(testis.list)<-"ident"
DimPlot(testis.list, split.by = "Age")


DimPlot(testis.list, label = T)
DimPlot(testis.list, label = T, split.by = "rep")

testis.list$Age<-ifelse(testis.list$rep %in% c("Old1","Old3a","Old3b"),"Old","Young")
Idents(testis.list)<-"Age"
allmarkers<-FindMarkers(testis.list, ident.1 = "Old")

testis.list$ident<-Idents(testis.list)
metadata<-testis.list@meta.data




#TE, denovo

Old1<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/Other data/210728_Old1_TE_denovo_run/")
Old1 <- CreateSeuratObject(Old1, min.cells = 3, min.features =200)
Old1$rep<-"Old1"


Old3a<-Read10X("Old3a_TE_denovo_run/")
Old3a <- CreateSeuratObject(Old3a, min.cells = 3, min.features =200)
Old3a$rep<-"Old3a"


Old3b<-Read10X("Old3b_TE_denovo_run/")
Old3b <- CreateSeuratObject(Old3b, min.cells = 3, min.features =200)
Old3b$rep<-"Old3b"




Rep1<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/Other data/210728_Young1_TE_denovo_run/")
Rep1 <- CreateSeuratObject(Rep1,min.cells = 3, min.features =200)
Rep1$rep<-"Rep1"



Rep2<-Read10X("~/Dropbox (Dropbox @RU)/Migrated From Box/Other data/210728_Young2_TE_denovo_run/")
Rep2 <- CreateSeuratObject(Rep2, min.cells = 3,min.features =200)
Rep2$rep<-"Rep2"


Younglib3<-Read10X("YoungLib3_TE_denovo_run/")
Younglib3 <- CreateSeuratObject(Younglib3, min.cells = 3, min.features =200)
Younglib3$rep<-"Younglib3"


testis.list2<-lapply(X = list(Old1, Old3a,Old3b, Rep1, Rep2, Younglib3), FUN = function(x){
  # x<-SCTransform(x)
  #x<-RunPCA(x)
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, nfeatures = 2000)
  x<-ScaleData(x)
  x<-RunPCA(x)
  # x<-RunUMAP(x, dims = 1:10)
  x<-SCTransform(x)
})


library(stringr)
features <- SelectIntegrationFeatures(object.list = testis.list2)
anchors <- FindIntegrationAnchors(object.list = testis.list2, anchor.features = features)
# this command creates an 'integrated' data assay
testis.list2 <- IntegrateData(anchorset = anchors)
#testis.list2<-readRDS("~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps_te_denovo.RDS")

testis.list2$barcode<-rownames(testis.list2@meta.data)

testis.list2$barcode<-stringr::str_split(string = testis.list2$barcode, pattern ="-", simplify = T)[,1]
metadata$barcode<-rownames(metadata)
metadata$barcode<-stringr::str_split(string = metadata$barcode, pattern ="-", simplify = T)[,1]
#join idents from other seurat object#####




metadata2<-join(testis.list2@meta.data, metadata, by = c("rep","barcode") )
testis.list2$ident<-metadata2$ident
#saveRDS(testis.list, "~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps_te_denovo.RDS")
testis.list2<-readRDS("~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220607_6reps_te_denovo.RDS")








####DE of de novo genes, TE #####

seg<-(scan("~/Dropbox (Dropbox @RU)/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/181221_li_segdenovo_DE_pseudotime.names.txt", what="", sep="\n"))
fixed<-(scan("~/Dropbox (Dropbox @RU)/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/181221_li_fixeddenovo_DE_pseudotime.names.txt", what="", sep="\n"))
seg <-unique(seg)
fixed <-unique(fixed)
denovo <- append(seg, fixed)
denovo <-unique(denovo)
denovo<-gsub(x=denovo, pattern ="_", replacement = "-" )
seg<-gsub(x=seg, pattern ="_", replacement = "-" )
fixed<-gsub(x=fixed, pattern ="_", replacement = "-" )

notdenovo<- testis.list2@assays[["RNA"]]@counts@Dimnames[[1]][!testis.list2@assays[["RNA"]]@counts@Dimnames[[1]] %in% denovo]


Idents(testis.list2)<- "ident"
avgexp<-AverageExpression(testis.list2)
avgexp<-avgexp$RNA



testis.list2<-ScaleData(testis.list2, assay  ="RNA", split.by = "Age")#testis.list2<-ScaleData(testis.list2, assay  ="RNA", split.by = "rep")
testis.list2@assays[["RNA"]]@scale.data[1:5,1:5]
testis.list2<-SCTransform(testis.list2)
testis.list2$Age<-ifelse(testis.list2$rep %in% c("Old1","Old3a","Old3b"),"Old","Young")


median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

library(ggpubr)

Idents(testis.list2)<-testis.list2$ident
tmp<-data.frame()
for ( i in c("Young", "Old")){
  for (j in na.omit(unique(testis.list2$ident))) {
    tmp2<-subset(testis.list2, Age == i & ident == j )@assays$RNA@scale.data
    tmp3<-data.frame(gene = rownames(tmp2), exp = rowMeans(tmp2), age = i, celltype = j)
    tmp<-rbind(tmp, tmp3)
  }
}

tmp$celltype<-factor(x = tmp$celltype, levels=c("Hub cells", "Cyst cells", "Epithelial cells","Accessory gland", "GSC, Early Spermatogonia", "Late spermatogonia", "Early spermatocytes", "Late spermatocytes", "Early spermatids", "Late spermatids"))
tmp<-subset(tmp, !celltype %in% c("Hub cells", "Epithelial cells", "Cyst cells", "Accessory gland"))
tmp<-na.omit(tmp)
tmp$Class<-ifelse(tmp$gene %in%  denovo, "de novo\n genes", "Other\ngenes")



TE<-read.delim("~/Dropbox (Dropbox @RU)/Migrated From Box/Other data/TE_list.txt", header = F)
TE<-TE$V1
tmp2<-subset(tmp, !gene %in% TE)
#keep
p1<-ggplot(tm2, aes(x= Class, y = exp, fill = age)) + ylim(-1,1)+
  theme_classic()+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~celltype, nrow = 1)+theme(axis.text=element_text(size=10))+stat_compare_means(method = "wilcox",p.adjust.methods = "bonferroni", label = "p.format")+scale_fill_manual(c("Dark"))




Idents(testis.list2)<-factor(Idents(testis.list2),levels=c("Hub cells", "Cyst cells", "Epithelial cells", "Accessory gland", "GSC, Early spermatogonia","Late spermatogonia", "Early spermatocytes", "Late spermatocytes", "Early spermatids", "Late spermatids"))
DimPlot(testis.list2, split.by = "Age")
#now do TEs


tmp<-subset(tmp, !gene %in% denovo)
tmp$Class<-ifelse(tmp$gene %in%  TE, "Transposable\n element", "Other\n genes")
tmp$Class<-factor(tmp$Class, levels = c("Transposable\n element","Other\n genes"))
tmp$age<-factor(tmp$age, levels = c("Young","Old"))
ggplot(tmp, aes(x= Class, y = exp, fill = age)) + ylim(-1,1)+
  theme_classic()+geom_boxplot(outlier.shape = NA)+ylab("scaled expression")+
  facet_wrap(~celltype, nrow = 1)+theme(axis.text=element_text(size=10))+stat_compare_means(method = "wilcox",p.adjust.methods = "bonferroni", label = "p.signif")
p2<-ggplot(tmp, aes(x= Class, y = exp, fill = age)) + ylim(-1,1)+
  theme_classic()+geom_boxplot(outlier.shape = NA)+ylab("scaled expression")+
  facet_wrap(~celltype, nrow = 1)+theme(axis.text=element_text(size=10))+xlab("")+stat_compare_means(method = "wilcox",p.adjust.methods = "bonferroni", label = "p.format")

p1<-p1+scale_fill_manual(values=c("royalblue","darkorchid4"))+ylab("scaled expression")+xlab("")
p2<-p2+scale_fill_manual(values=c("royalblue","darkorchid4"))
#figure 5
plot_grid(p1,p2, labels = c("A","B"),nrow = 2)
ggsave("Figure5.pdf", width= 14, height = 8)
#message- post-meiotic upregulation



#TE markers, de novo markers

Idents(testis.list2)<-"Age"
markers<-FindAllMarkers(testis.list2, only.pos = T, assay  = "RNA")


denovomarkers<-subset(data.frame(markers), gene %in% denovo)
TEmarkers<-subset(data.frame(markers), gene %in% TE)

# split between cell type and age


Idents(testis.list2)<-c("Age")
tmp<-data.frame()
for (i in na.omit(unique(testis.list2$ident))){
  tmp2<-FindAllMarkers(subset(testis.list2, ident == i), only.pos = T, assay = "SCT")
  tmp2$celltype<-i
  tmp<-rbind(tmp, tmp2)
}

denovomarkers<-subset(data.frame(tmp), gene %in% denovo)
segmarkers<-subset(denovomarkers, gene %in% seg)
segmarkers$type<-"seg"
fixedmarkers<-subset(denovomarkers, gene %in% fixed)
fixedmarkers$type<-"fixed"
denovomarkers<-rbind(segmarkers, fixedmarkers)
TEmarkers<-subset(data.frame(tmp), gene %in% TE)
#supplemental tables)
#write.csv(denovomarkers,"denovo.csv")
#write.csv(TEmarkers,"TE.csv")

#make plot with number of de novo, TE enriched in each cell type, young or old


tmp2<-unique(TEmarkers[,c("cluster","gene")])$cluster %>%table() %>% data.frame()
names(tmp2)<-c("cluster","Freq")
tmp2$type<-"Transposable Element"
tmp3<-denovomarkers %>%group_by(cluster, type) %>%tally()

tmp3<-unique(denovomarkers[,c("cluster","type","gene")])[,c("cluster","type")] %>%table() %>%data.frame()


tmp3<-rbind(tmp2, tmp3)


tmp3$type<-gsub(x=tmp3$type,  pattern ="seg", "Segregating de novo")
tmp3$type<-gsub(x=tmp3$type,pattern ="fixed", "Fixed de novo")
tmp3$cluster<-factor(tmp3$cluster, levels = c("Young","Old"))

#figure 5####
ggplot(tmp3, aes(x=cluster,fill= cluster, y=Freq))+geom_bar(stat= "identity")+facet_wrap(~type)+theme_classic()+scale_fill_manual(values=c("royalblue","darkorchid4"))+
  theme(legend.position = "none")+ylab("Number of elements enriched")+scale_y_continuous(breaks=c(0,2,4,6,8,10))
#ggsave("~/Dropbox (Dropbox @RU)/Evan/Zhao lab notebook/Misc_data/210625_figure5_seg_fixed_TE.pdf", width=7, height=3)
#reviewer didn't like it. let's do supp 3 as figure 5




#DN/DS by age, cell type


#####figure 6#######
testis.list$Age<-ifelse(testis.list$rep %in% c("Old1","Old3a","Old3b"),"Old","Young")


dnds<-read.delim("~/Downloads/melsubgroup_analysis_results_flydivas_v1.2")

dnds<-dnds[,c("symbol", "dN", "dS")]
colnames(dnds)<-c("gene","dN", "dS")

Idents(testis.list)<-"Age"
#get enriched  marker genes between cell types, ages
Agemarkers<-FindAllMarkers(testis.list, only.pos = T )
Agemarkers<-plyr::join(Agemarkers, dnds)
Agemarkers<-na.omit(Agemarkers)
Agemarkers$dnds<-Agemarkers$dN/Agemarkers$dS
ggplot(Agemarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.2)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("All cells")

gscmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="GSC, Early spermatogonia"),only.pos = T, assay = "RNA")
#write.csv(gscmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_GSC_EarlySpermatogonia_markers.csv")

latespgmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatogonia"),only.pos = T, assay = "RNA")
#write.csv(latespgmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_LateSpermatogonia_markers.csv")

earlyspermatocytemarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Early spermatocytes"),only.pos = T, assay = "RNA")
#write.csv(earlyspermatocytemarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Earlyspermatocyte_markers.csv")

latespermatocytemarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatocytes"),only.pos = T, assay = "RNA")
#write.csv(latespermatocytemarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Latespermatocyte_markers.csv")

earlyspermatidmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Early spermatids"),only.pos = T, assay = "RNA")
#write.csv(earlyspermatidmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Earlyspermatid_markers.csv")

latespermatidmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatids"),only.pos = T, assay = "RNA")
#write.csv(latespermatidmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Latespermatid_markers.csv")



gscmarkers<-plyr::join(gscmarkers, dnds)
gscmarkers<-na.omit(gscmarkers)
gscmarkers$dnds<-gscmarkers$dN/gscmarkers$dS
ggplot(gscmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("GSC, Early spermatogonia")

latespgmarkers<-plyr::join(latespgmarkers, dnds)
latespgmarkers<-na.omit(latespgmarkers)
latespgmarkers$dnds<-latespgmarkers$dN/latespgmarkers$dS
ggplot(latespgmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")



latespermatidmarkers<-plyr::join(latespermatidmarkers, dnds)
latespermatidmarkers<-na.omit(latespermatidmarkers)
latespermatidmarkers$dnds<-latespermatidmarkers$dN/latespermatidmarkers$dS
ggplot(latespermatidmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("Late spermatids")



gscmarkers$Celltype<-"GSC, Early spermatogonia"
latespermatidmarkers$Celltype<-"Late spermatids"
youngmarkers<-rbind(subset(gscmarkers, cluster =="Young"), subset(latespermatidmarkers, cluster == "Young"))
oldmarkers<-rbind(subset(gscmarkers, cluster =="Old"), subset(latespermatidmarkers, cluster == "Old"))

combinedmarkers<-rbind(oldmarkers, youngmarkers)

ggplot(combinedmarkers, aes(x=Celltype, y = dnds, fill = Celltype))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+facet_wrap(~cluster)+
  stat_summary(aes( x=Celltype, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")




latespermatocytemarkers<-plyr::join(latespermatocytemarkers, dnds)
latespermatocytemarkers<-na.omit(latespermatocytemarkers)
latespermatocytemarkers$dnds<-latespermatocytemarkers$dN/latespermatocytemarkers$dS
ggplot(latespermatocytemarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("Late spermatocytes")

Idents(testis.list)<-"ident"
earlymarkers1<-FindMarkers(subset(testis.list, age == "Young"), ident.1 = "GSC, Early spermatogonia",ident.2 = "Late spermatids", only.pos = T)
latemarkers1<-FindMarkers(subset(testis.list, age == "Young"), ident.2 = "GSC, Early spermatogonia",ident.1 = "Late spermatids", only.pos = T)
earlymarkers1$celltype<-"GSC, Early spermatogonia"
latemarkers1$celltype<-"Late spermatids"
earlymarkers1$age<-"Young"
latemarkers1$age<-"Young"
earlymarkers2<-FindMarkers(subset(testis.list, age == "Old"), ident.1 = "GSC, Early spermatogonia",ident.2 = "Late spermatids", only.pos = T)
latemarkers2<-FindMarkers(subset(testis.list, age == "Old"), ident.2 = "GSC, Early spermatogonia",ident.1 = "Late spermatids", only.pos = T)
earlymarkers2$celltype<-"GSC, Early spermatogonia"
latemarkers2$celltype<-"Late spermatids"
earlymarkers2$age<-"Old"
latemarkers2$age<-"Old"
earlymarkers1$gene<-rownames(earlymarkers1)
earlymarkers2$gene<-rownames(earlymarkers2)
latemarkers1$gene<-rownames(latemarkers1)
latemarkers2$gene<-rownames(latemarkers2)



earlylatemarkers<-do.call("rbind",list(earlymarkers1, latemarkers1, earlymarkers2, latemarkers2) )
earlylatemarkers$gene<-rownames(earlylatemarkers)
earlylatemarkers<-plyr::join(earlylatemarkers,dnds)
earlylatemarkers<-na.omit(earlylatemarkers)
earlylatemarkers$dnds<-earlylatemarkers$dN/earlylatemarkers$dS
earlylatemarkers$age<-factor(earlylatemarkers$age, levels = c("Young","Old"))

#get p values
stat.test <- compare_means(
  dnds ~ celltype, data = earlylatemarkers, group.by = "age",
  method = "wilcox"
)


#figure 6A
p1<-ggplot(subset(earlylatemarkers,celltype %in% c("GSC, Early spermatogonia","Late spermatids")), aes(x=celltype, y = dnds, fill = age))+facet_wrap(~age)+geom_violin()+theme_classic()+
  stat_summary(aes( x=celltype, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("")+
  theme(axis.text=element_text(size = 10))+
  ggtitle("Cell type-enriched genes")+
  theme(title=element_text(size=15),axis.text=element_text(size = 15), axis.title = element_text(size=15), strip.text = element_text(size=15), legend.position = "none")+ 
  stat_pvalue_manual(data = stat.test,bracket.shorten = 1,tip.length = 0, y.position = 0.6)+scale_fill_manual(values=c("royalblue","darkorchid4"))

#getr



Idents(testis.list)<-"age"
#get enriched  marker genes between cell types, ages
Agemarkers<-FindAllMarkers(testis.list, only.pos = T )
#write.csv(Agemarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Agemarkers.csv")
Idents(testis.list)<-"age"

gscmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="GSC, Early spermatogonia"),only.pos = T, assay = "RNA")
#write.csv(gscmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_GSC_EarlySpermatogonia_markers.csv")

latespgmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatogonia"),only.pos = T, assay = "RNA")
#write.csv(latespgmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_LateSpermatogonia_markers.csv")

earlyspermatocytemarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Early spermatocytes"),only.pos = T, assay = "RNA")
#write.csv(earlyspermatocytemarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Earlyspermatocyte_markers.csv")

latespermatocytemarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatocytes"),only.pos = T, assay = "RNA")
#write.csv(latespermatocytemarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Latespermatocyte_markers.csv")

earlyspermatidmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Early spermatids"),only.pos = T, assay = "RNA")
#write.csv(earlyspermatidmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Earlyspermatid_markers.csv")

latespermatidmarkers<-FindAllMarkers(subset(testis.list, subset = ident =="Late spermatids"),only.pos = T, assay = "RNA")
#write.csv(latespermatidmarkers, "~/Box Sync/Evan/Zhao lab notebook/Misc_data/210603_R517_Latespermatid_markers.csv")

gscmarkers<-plyr::join(gscmarkers, dnds)
gscmarkers<-na.omit(gscmarkers)
gscmarkers$dnds<-gscmarkers$dN/gscmarkers$dS
ggplot(gscmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("GSC, Early spermatogonia")

latespgmarkers<-plyr::join(latespgmarkers, dnds)
latespgmarkers<-na.omit(latespgmarkers)
latespgmarkers$dnds<-latespgmarkers$dN/latespgmarkers$dS
ggplot(latespgmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")



latespermatidmarkers<-plyr::join(latespermatidmarkers, dnds)
latespermatidmarkers<-na.omit(latespermatidmarkers)
latespermatidmarkers$dnds<-latespermatidmarkers$dN/latespermatidmarkers$dS
ggplot(latespermatidmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(label.x = 1.5)+theme_classic()+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("Fly age")+scale_fill_manual(values=c("royalblue","grey"))+theme(axis.text=element_text(size = 10), legend.position = "none")+ggtitle("Late spermatids")



gscmarkers$Celltype<-"GSC, Early spermatogonia"
latespermatidmarkers$Celltype<-"Late spermatids"
youngmarkers<-rbind(subset(gscmarkers, cluster =="Young"), subset(latespermatidmarkers, cluster == "Young"))
oldmarkers<-rbind(subset(gscmarkers, cluster =="Old"), subset(latespermatidmarkers, cluster == "Old"))

combinedmarkers<-rbind(oldmarkers, youngmarkers)

ggplot(combinedmarkers, aes(x=Celltype, y = dnds, fill = cluster))+geom_violin()+stat_compare_means(aes(fill = Celltype),label.x = 1.5)+theme_classic()+facet_wrap(~cluster)+
  stat_summary(aes( x=Celltype, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+scale_fill_manual(values=c("royalblue","darkorchid4"))+theme(axis.text=element_text(size = 15), axis.title = element_text(size=15), strip.text = element_text(size=15), legend.position = "none")

stat.test2 <- combinedmarkers %>%
  group_by(Celltype) %>%
  wilcox_test(dnds ~ cluster) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>% add_y_position(fun ="max")


stat.test2 <- combinedmarkers %>%
  group_by(Celltype) %>%
  wilcox_test(dnds ~ cluster) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test2$y.position[]<-1
stat.test2$p.adj<-formatC(stat.test2$p.adj, digits = 2, format = "e")
combinedmarkers$cluster<-factor(combinedmarkers$cluster, levels=c("Young","Old"))


p2<-ggplot(combinedmarkers, aes(x=cluster, y = dnds, fill = cluster))+geom_violin()+theme_classic()+facet_wrap(~Celltype)+
  stat_summary(aes( x=cluster, y=dnds),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5
  )+ylab("dN/dS")+xlab("")+scale_fill_manual(values=c("royalblue","darkorchid4"))+
  ggtitle("Age-enriched genes")+
  theme(title=element_text(size=15),
        axis.text=element_text(size = 15), 
        axis.title = element_text(size=15), strip.text = element_text(size=15),
        legend.position = "none")+
            geom_text( aes(x,y,label = lab), 
                       data = data.frame(x=c("Young","Young"), y=0.6, cluster = rep(c("Old","Young"),2),
                                         Celltype = stat.test2$Celltype,lab=stat.test2$p.adj), hjust = -0.45)





#Figure 6, now figure 5

plot_grid(p1, p2, nrow = 2, labels = c("A","B"), label_size = 15)
ggsave("Figure_5.pdf", width=10, height=8)


#figure 3######

testis.list$Age<-ifelse(testis.list$rep %in% c("Old1","Old3a","Old3b"),"Old","Young")
testis.list$ident<-Idents(testis.list)
DDR<- (scan("~/Dropbox (Dropbox @RU)/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/DDR2.txt", what="", sep="\n" ))
testis.list<-ScaleData(testis.list, assay  ="SCT", split.by = "Age")
testis.list@assays[["RNA"]]@scale.data[1:5,1:5]


median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}


tmp<-data.frame()
for ( i in c("Old", "Young")){
  for (j in unique(testis.list$ident)) {
    tmp2<-subset(testis.list, Age == i & ident == j )@assays$RNA@scale.data
    tmp3<-data.frame(gene = rownames(tmp2), exp = rowMeans(tmp2), Age = i, celltype = j)
    tmp<-rbind(tmp, tmp3)
  }
}

tmp$celltype<-factor(x = tmp$celltype, levels=c("Hub cells", "Cyst cells", "Epithelial cells", "GSC, Early spermatogonia", "Accessory gland", "Late spermatogonia", "Early spermatocytes", "Late spermatocytes", "Early spermatids", "Late spermatids"))
tmp<-subset(tmp, !celltype %in% c("Hub cells", "Epithelial cells", "Cyst cells"))
tmp$DDR<-ifelse(tmp$gene %in%  DDR, "Genome\nmaintenance\ngenes", "Other\ngenes")
ggplot(tmp, aes(x= celltype, y = exp, fill = DDR)) + ylim(-1,1)+
  theme_classic()+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Age, nrow = 2)+theme(axis.text=element_text(size=10))



#add p values
library(rstatix)
stat.test <- na.omit(tmp) %>%
  group_by(celltype,DDR) %>%
  wilcox_test(exp ~ Age) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>% add_y_position(fun ="max")
stat.test$y.position[]<-1
stat.test$p.adj<-formatC(stat.test$p.adj, digits = 2, format = "e")
tmp$Age<-factor(tmp$Age, levels=c("Young","Old"))
#best one
fig3b<-ggplot(na.omit(tmp), aes(x= DDR, y = exp)) + ylim(-1,1)+
  theme_classic()+geom_violin(scale = "width",aes(fill = Age))+labs(fill = "Fly age", y = "Scaled expression",x= "Gene type")+
  facet_wrap(~celltype, nrow = 1)+theme(axis.text=element_text(size=10))+
  stat_summary(aes( x=DDR, y=exp, fill = Age),fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5 )+
  geom_text(aes(x, y, label=lab), data=data.frame(x = rep(c("Genome\nmaintenance\ngenes","Other\ngenes"),6), y=Inf, celltype = stat.test$celltype, lab=stat.test$p.adj, yy=letters[1:3]), vjust=1)+
  scale_fill_manual(values=c("royalblue","darkorchid4"))+xlab("")



#Figure 4 #####

options(ggrepel.max.overlaps=Inf)
gsc <-subset(testis.list, subset = ident == "GSC, Early spermatogonia")
Idents(gsc) <- "Age"
late.spermatogonia <-subset(testis.list, subset = ident == "Late spermatogonia")
Idents(late.spermatogonia) <- "Age"
early.spermatocytes <-subset(testis.list, subset = ident == "Early spermatocytes")
Idents(early.spermatocytes)<-"Age"
late.spermatocytes <-subset(testis.list, subset = ident == "Late spermatocytes")
Idents(late.spermatocytes)<-"Age"
early.spermatids <-subset(testis.list, subset = ident == "Early spermatids")
Idents(early.spermatids)<-"Age"
late.spermatids<-subset(testis.list, subset = ident == "Late spermatids")
Idents(late.spermatids)<-"Age"

gsc.age.markers<-FindMarkers(gsc, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
gsc.age.markers$gene <-rownames(gsc.age.markers)
genes.label<-subset(gsc.age.markers, p_val_adj <.05  & gene %in% DDR )$gene
keyvals <- ifelse(
  gsc.age.markers$avg_log2FC < 0 & gsc.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(gsc.age.markers$avg_log2FC > 0 & gsc.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'

p1<-EnhancedVolcano(gsc.age.markers,labSize = 3, title = "GSC, Early spermatogonia", drawConnectors = T,
                    lab = gsc.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)

p1

late.spermatogonia.age.markers<-FindMarkers(late.spermatogonia, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
late.spermatogonia.age.markers$gene <-rownames(late.spermatogonia.age.markers)
genes.label<-subset(late.spermatogonia.age.markers, p_val_adj <.05  & gene %in% DDR )$gene
keyvals <- ifelse(
  late.spermatogonia.age.markers$avg_log2FC < 0 & late.spermatogonia.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(late.spermatogonia.age.markers$avg_log2FC > 0 & late.spermatogonia.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'
p2<-EnhancedVolcano(late.spermatogonia.age.markers,labSize = 3, title = "Late spermatogonia", drawConnectors = T,
                    lab = late.spermatogonia.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)


  
  
  early.spermatocytes <-subset(testis.list, subset = ident == "Early spermatocytes")
Idents(early.spermatocytes) <- "Age"

early.spermatocytes.age.markers<-FindMarkers(early.spermatocytes, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
early.spermatocytes.age.markers$gene <-rownames(early.spermatocytes.age.markers)
genes.label<-subset(early.spermatocytes.age.markers, p_val_adj <.05 & abs(avg_log2FC) >1 & gene %in% DDR )$gene
genes.label<-subset(late.spermatogonia.age.markers, p_val_adj <.05  & gene %in% DDR )$gene
keyvals <- ifelse(
  early.spermatocytes.age.markers$avg_log2FC < 0 & early.spermatocytes.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(early.spermatocytes.age.markers$avg_log2FC > 0 & early.spermatocytes.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'

p3<-EnhancedVolcano(early.spermatocytes.age.markers,labSize = 3, title = "Early spermatocytes", drawConnectors = T,
                    lab = early.spermatocytes.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)
p3


late.spermatocytes <-subset(testis.list, subset = ident == "Late spermatocytes")
Idents(late.spermatocytes) <- "Age"

late.spermatocytes.age.markers<-FindMarkers(late.spermatocytes, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
late.spermatocytes.age.markers$gene <-rownames(late.spermatocytes.age.markers)

genes.label<-subset(late.spermatocytes.age.markers, p_val_adj <.05 & abs(avg_log2FC) >1 & gene %in% DDR )$gene
keyvals <- ifelse(
  late.spermatocytes.age.markers$avg_log2FC < 0 & late.spermatocytes.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(late.spermatocytes.age.markers$avg_log2FC > 0 & late.spermatocytes.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'

p4<-EnhancedVolcano(late.spermatocytes.age.markers,labSize = 3, title = "Late spermatocytes", drawConnectors = T,
                    lab = late.spermatocytes.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)
p4



early.spermatids <-subset(testis.list, subset = ident == "Early spermatids")
Idents(early.spermatids) <- "Age"

early.spermatids.age.markers<-FindMarkers(early.spermatids, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
early.spermatids.age.markers$gene <-rownames(early.spermatids.age.markers)

genes.label<-subset(early.spermatids.age.markers, p_val_adj <.05 & abs(avg_log2FC) >1 & gene %in% DDR )$gene
keyvals <- ifelse(
  early.spermatids.age.markers$avg_log2FC < 0 & early.spermatids.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(early.spermatids.age.markers$avg_log2FC > 0 & early.spermatids.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'

p5<-EnhancedVolcano(early.spermatids.age.markers,labSize = 3, title = "Early spermatids", drawConnectors = T,
                    lab = early.spermatids.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)
p5


late.spermatids <-subset(testis.list, subset = ident == "Late spermatids")
Idents(late.spermatids) <- "Age"

late.spermatids.age.markers<-FindMarkers(late.spermatids, ident.1 = "Young", ident.2 = "Old", assay = "SCT")
late.spermatids.age.markers$gene <-rownames(late.spermatids.age.markers)

genes.label<-subset(late.spermatids.age.markers, p_val_adj <.05 & abs(avg_log2FC) >1 & gene %in% DDR )$gene
keyvals <- ifelse(
  late.spermatids.age.markers$avg_log2FC < 0 & late.spermatids.age.markers$p_val_adj <.05, 'darkorchid4',
  ifelse(late.spermatids.age.markers$avg_log2FC > 0 & late.spermatids.age.markers$p_val_adj<.05, 'royalblue',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Young enriched'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == 'darkorchid4'] <- 'Old enriched'

p6<-EnhancedVolcano(late.spermatids.age.markers,labSize = 3, title = "Late spermatids", drawConnectors = T,
                    lab = late.spermatids.age.markers$gene, selectLab= genes.label,colConnectors = "red",endsConnectors = "last",
                    x = 'avg_log2FC', pCutoff = 0.05, FCcutoff = 0,colCustom = keyvals,boxedLabels = F,subtitle  = NULL,colAlpha = 0.4,
                    y = 'p_val_adj', caption = NULL)+xlim(-2.5,2.5)

p6


#flip x axis
p1<-p1+scale_x_reverse()
p2<-p2+scale_x_reverse()
p3<-p3+scale_x_reverse()
p4<-p4+scale_x_reverse()
p5<-p5+scale_x_reverse()
p6<-p6+scale_x_reverse()+theme(legend.position = 'none')

p7<-plot_grid(p1,p2,p3,p4,p5,p6, nrow = 3, labels = c("A","B","C","D","E","F"), label_size = 20)
ggsave(plot = p7,"Figure_4.pdf", width=20, height=25)
#table 1: numbers of genes

for (i in ls()[grep("age.markers",ls())]){
  print(i)
  print(nrow(subset(get(i), avg_log2FC>0)))
  print(nrow(subset(get(i), avg_log2FC<0)))
  
}

#supplemental table 4: enrichment statistics of genes by cell type
tmp<-data.frame()
for (i in ls()[grep("age.markers",ls())]){
  print(i)
  tmp2<-data.frame(subset(get(i), p_val_adj<.05 & gene %in% DDR) )
  #tmp2<-data.frame(subset(get(i), p_val_adj<.05) )
  tmp2$celltype<-i
  tmp2$celltype<-gsub(x=tmp2$celltype,pattern = "age.markers","")
  tmp2$celltype<-gsub(x=tmp2$celltype,pattern = '\\.'," ")
  tmp<-rbind(tmp, tmp2)
  
}

tmp$Age<-ifelse(tmp$avg_log2FC>0, "Young","Old")
write.csv(tmp, "All_markers.csv")

#supplemental table 1#####


tmp2<-data.frame()
tmp3<-data.frame(subset(gsc.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"GSC, Early Spermtatogonia"
tmp2<-rbind(tmp2, tmp3)

tmp3<-data.frame(subset(late.spermatogonia.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"Late spermatogonia"
tmp2<-rbind(tmp2, tmp3)

tmp3<-data.frame(subset(early.spermatocytes.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"Early spermatocytes"
tmp2<-rbind(tmp2, tmp3)

tmp3<-data.frame(subset(late.spermatocytes.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"Late spermatocytes"
tmp2<-rbind(tmp2, tmp3)

tmp3<-data.frame(subset(early.spermatids.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"Early spermatids"
tmp2<-rbind(tmp2, tmp3)

tmp3<-data.frame(subset(late.spermatids.age.markers, p_val <.05 & gene %in% DDR))
tmp3$celltype <-"Late spermatids"
tmp2<-rbind(tmp2, tmp3)

tmp2$Age<-ifelse(tmp2$avg_log2FC>0, "Young","Old")

####panther####
Idents(testis.list)<-"Age"
allmarkers<-FindAllMarkers(testis.list, assay = "SCT")

noquote(subset(allmarkers, p_val < .05 & avg_log2FC > 0)$gene)
noquote(subset(allmarkers, p_val < .05 & avg_log2FC < 0)$gene)
library(stringr)

###Figure 2
library(stringr)
se <- function(x) sqrt(var(x) / length(x))
SNPS<-read.csv("~/Dropbox (Dropbox @RU)/220605_Mutation_Paper_10X_Cellranger/220628_SNPS_readnames_barcodes_coverage.csv")


SNPS2<-read.csv("~/Dropbox (Dropbox @RU)/Migrated From Box/R517_3reps/210607_SNPS_readnames_barcodes_coverage.csv")
SNPS$barcode<-str_split(SNPS$barcode, pattern = ":", simplify = T)[,3]
SNPS$barcode<-str_split(SNPS$barcode, pattern = "-", simplify = T)[,1]
SNPS2$barcode<-str_split(SNPS2$barcode, pattern = ":", simplify = T)[,3]
SNPS2$barcode<-str_split(SNPS2$barcode, pattern = "-", simplify = T)[,1]
SNPS<-rbind(SNPS, SNPS2)
testis.list$barcode<-rownames(testis.list@meta.data)
testis.list$barcode<-str_split(string = testis.list$barcode, pattern ="-", simplify = T)[,1]
SNPS<-join(SNPS, testis.list@meta.data, by = c("rep","barcode"))

SNPS<-subset(SNPS, !ident =="NA") #this is fixed now
SNPS<-subset(SNPS, coverage > 10) 
head(SNPS)
tail(SNPS)


somatic <-subset(SNPS, ident %in% c("Cyst cells", "Hub cells", "Epithelial cells", "Accessory gland"))
somatic <-unique(somatic[,c("pos","chr", "barcode","rep")])[,c("pos","chr")]
somatic<-paste(somatic$pos, somatic$chr, somatic$rep,sep = ",")
somatic<-somatic %>%table() %>%data.frame() %>%subset(Freq > 0)
SNPS$.<-paste(SNPS$pos, SNPS$chr, SNPS$rep,sep = ",")


 
SNPS<-subset(SNPS, !. %in% somatic$.)
#

#keep SNPS present in only one rep, count UMIs per cell
SNPS$mutation<-paste(SNPS$ref, SNPS$alt)
head(SNPS)
data5<-SNPS %>% mutate (merged = paste(chr, pos,ref,alt)) 
head(data5)
data5<-unique(data5[,c("merged","rep")])[,c("merged")] %>%table() %>%data.frame() %>%subset(Freq ==1)
head(data5)


#library(stringr)
data5<-stringr::str_split(data5$., pattern=" ", simplify=T)%>%data.frame()
head(data5)
colnames(data5)<-c("chr", "pos","ref","alt")
data5<-join(data5, SNPS)
head(data5)


unique(data5[,c("chr", "pos", "ref", "alt","mutation","Age")])[,c("Age","mutation")] %>%table() %>%data.frame()
#this gets the proportions used for the next part
unique(data5[,c("chr", "pos", "ref", "alt","mutation","Age")])[,c("Age","mutation")] %>%table() %>%data.frame() %>% group_by(Age) %>%summarise(sum=sum(Freq))


#merge clustered mutations
data3<-unique(data5[,c("chr","pos", "ref", "alt")])
data3$pos<-as.numeric(data3$pos)
data3<-data3 %>% group_by(chr) %>%mutate(pos_bin = cut(pos, breaks = seq(0,max(pos), by = 10)))
data3$pos_bin<-as.character(data3$pos_bin)

data3<-data3 %>% 
  dplyr::group_by(chr, pos_bin) %>% 
  dplyr::slice(which.min(pos))

data5<-join(data3, data5)

#ensure that no SNPS are present in multiple reps:
 nrow(unique(data5[,c("chr","pos","ref","alt", "rep")]))
#5336
nrow(unique(data5[,c("chr","pos","ref","alt")]))
#5336
#same number, means no dups.

#require more than one read:
data5<-data5 %>% group_by(chr,pos, ref, alt,  rep) %>% mutate (nMutReads= n_distinct(readname))
data5<-subset(data5, nMutReads > 1)

data5$Age<-ifelse(data5$rep %in% c("Old1","Old3a","Old3b"),"Old", "Young")
#get correct number of cells. get number of unique barcodes matched to mutations. divide from n

tmp1<- unique(data5[,c("barcode","ident","rep")])[,c("barcode","ident","rep")] %>% group_by(rep, ident) %>% mutate(nMutCell=n_distinct(barcode)) %>% data.frame()


tmp2<-testis.list@meta.data %>% group_by(rep, ident) %>% mutate(n=n_distinct(barcode), meanreads = mean(nCount_RNA)) %>%data.frame()
tmp3<-join(tmp1, tmp2)
tmp3$prop<-tmp3$nMutCell/tmp3$n
#tmp3$merged<-tmp3$nMutCell/(tmp3$n*tmp3$meanreads)


tmp3$Age<-ifelse(tmp3$rep %in% c("Old1","Old3a","Old3b"),"Old", "Young")
tmp3<-subset(tmp3, ident %in% c("GSC, Early spermatogonia","Late spermatogonia","Early spermatocytes","Late spermatocytes",
                                "Early spermatids","Late spermatids"))

tmp3$ident<-factor(tmp3$ident, levels = c("GSC, Early spermatogonia","Late spermatogonia","Early spermatocytes","Late spermatocytes",
                                          "Early spermatids","Late spermatids"))

#alt figure 2: >10 ref reads somatic not allowed,need to add error bars
tmp3<-tmp3 %>%group_by(Age, ident) %>%dplyr::mutate(avgprop=mean(unique(prop)))
tmp3<-tmp3 %>% group_by(Age, ident, avgprop) %>% dplyr::mutate(ymin=avgprop-se(unique(prop)), ymax = avgprop+se(unique(prop)))

tmp3$Age<-factor(tmp3$Age, levels=c("Young","Old"))

#add p values
tmp4<-data.frame()
  for (i in levels(tmp3$ident)){
    tmp5<-subset(tmp3,  ident == i)
    tmp5<-unique(tmp5[,c("Age","rep","prop","n")])
    tmp6<-subset(tmp5, Age =="Young") %>%summarise(meanprop=mean(prop), n=sum(n))
    tmp7<-subset(tmp5, Age =="Old") %>%summarise(meanprop=mean(prop), n=sum(n))
    p<-prop.test(c(tmp6$meanprop*tmp6$n, tmp7$meanprop*tmp7$n), c(tmp6$n, tmp7$n))$p.value
    #p<-wilcox.test(subset(tmp5, Age =="Young")$prop,subset(tmp5, Age =="Old")$prop  )$p.value
    tmp4<-rbind(tmp4, data.frame(ident = i, p = p))
    }


tmp4$p.adjust<-p.adjust(tmp4$p)
tmp4$p.adjust<-formatC(tmp4$p.adjust, digits = 2, format = "e")
tmp3<-join(tmp3, tmp4)



fig2alt<-ggplot(unique(tmp3[,c("ident","avgprop","Age","ymin","ymax","p.adjust")]), aes(x = ident, y = avgprop, col = Age, fill = Age))+geom_bar(position = "dodge",stat = "identity")+
  theme_classic()+theme(axis.text=element_text(size=10) )+
  geom_errorbar(width=0.3, position= position_dodge(width=0.9),col= "black",size=0.5,aes(ymin=ymin,ymax=ymax))+
  labs(fill = "Fly age", col = "Fly age")+
  scale_fill_manual(values=c("royalblue","darkorchid4"))+
  scale_color_manual(values=c("royalblue","darkorchid4"))+
  ylab("Proportion of cells with mutations")+xlab("Cell type")+geom_text(y=1,aes(x=ident,label = p.adjust), col = "black",check_overlap = T )





#count average number of SNPS in a given cell type- try splitting by cell?
tmp1<- unique(data5[,c("rep","ident","ref","alt","pos","barcode","nCount_RNA")]) %>% group_by(rep, ident, barcode) %>% summarise(nMut=n()) %>% data.frame()

#calculate proportions individually by rep
#tmp2<-testis.list@meta.data %>% group_by(rep, ident,barcode) %>% dplyr::summarise(n=n_distinct(barcode), meanreads = mean(nCount_RNA)) %>%data.frame()
tmp2<-testis.list@meta.data
tmp3<-join(tmp1, tmp2, by = c("rep","ident","barcode"))
tmp3$prop<-tmp3$nMut/tmp3$nCount_RNA

#tmp3$merged<-tmp3$Freq/(tmp3$n*tmp3$sumreads)


tmp3$Age<-ifelse(tmp3$rep %in% c("Old1","Old3a","Old3b"),"Old", "Young")
tmp3<-subset(tmp3, ident %in% c("GSC, Early spermatogonia","Late spermatogonia","Early spermatocytes","Late spermatocytes",
                                "Early spermatids","Late spermatids"))

tmp3$ident<-factor(tmp3$ident, levels = c("GSC, Early spermatogonia","Late spermatogonia","Early spermatocytes","Late spermatocytes",
                                          "Early spermatids","Late spermatids"))
#tmp3<-tmp3 %>% group_by(Age, ident) %>% mutate(avgmut=mean(merged))


tmp3<-tmp3 %>% group_by(Age, ident) %>% dplyr::mutate(SE=se(prop))
tmp3<-tmp3 %>% group_by(Age, ident) %>% dplyr::mutate(meanprop=mean(prop))

tmp3$Age<-factor(tmp3$Age, levels=c("Young","Old"))




#do wilcox test on mean snpsperumi. check p values.
tmp5<-data.frame()
for (i in levels(tmp3$ident)){
  p<-wilcox.test(unique(subset(tmp3, ident ==i & Age =="Young")$prop), unique(subset(tmp3, ident ==i & Age =="Old")$prop))$p.value
  p<-formatC(p, digits = 2, format = "e")
  tmp5<-rbind(tmp5, data.frame(p=p,ident = i))
}
tmp5$p.adjust<-p.adjust(tmp5$p)
tmp5$p.adjust<-formatC(tmp5$p.adjust, digits = 2, format = "e")
tmp3<-join(tmp3, tmp5, by = "ident")
tmp3$p.adjust

#try wilcox
wilcox.test(subset(tmp3, ident =="Late spermatids" & Age =="Young")$prop, subset(tmp3, ident =="Late spermatids" & Age =="Old")$prop)


p2<-ggplot(tmp3, aes(x = ident, y = prop))+geom_violin(scale = "width",aes(fill =Age ))+
  theme_classic()+theme(axis.text=element_text(size=10) )+  labs(fill = "Fly age", col = "Fly age")+scale_fill_manual(values=c("royalblue","darkorchid4"))+scale_color_manual(values=c("royalblue","darkorchid4"))+
  ylab("SNPS per UMI")+xlab("Cell type")+ylim(0,.005)+
  geom_text(check_overlap = T,y = .004,aes(x = ident,  label = p.adjust))+
  stat_summary(aes( x=ident, y=prop, fill = Age),
               fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5 )



plot_grid(fig2alt, p2, nrow = 2)
ggsave("Figure2.pdf", width= 12, height=8)




write.csv(unique(data5), "210603_R517_Old_Young_SNP_database.csv", col.names = T)



#make summary table of nCell, numi

tmp<-testis.list@meta.data %>% group_by(rep)%>%summarise(nCell= n(), umipercell=mean(nCount_RNA))


#figure 3: mutational signatures######



#total mutations in cell type
#for each cell type, count all the unique mutations of each type. divide each by number of total different mutations detected in that cell type. 
data6<-unique(data5[,c("mutation","Age","barcode","ident", "rep")]) 
data6<-na.omit(data6)
data6$mutation<-gsub(x=data6$mutation,pattern = "A C",replacement = "T G")
data6$mutation<-gsub(x=data6$mutation,pattern = "A G",replacement = "T C")
data6$mutation<-gsub(x=data6$mutation,pattern = "A T",replacement = "T A")
data6$mutation<-gsub(x=data6$mutation,pattern = "G T",replacement = "C A")
data6$mutation<-gsub(x=data6$mutation,pattern = "G C",replacement = "C G")
data6$mutation<-gsub(x=data6$mutation,pattern = "G A",replacement = "C T")


#get p values- each set of props
data6<-table(data6[,c("mutation", "rep","Age")])%>%data.frame() %>%group_by(rep) %>%mutate(prop=Freq/sum(Freq))
data6<-subset(data6, Freq > 0)
data6<-data6 %>%group_by(Age, mutation) %>%dplyr::mutate(meanprop=mean(prop), seprop=se(prop))
data6<-data6 %>%group_by(Age, mutation) %>%mutate(meanfreq=mean(Freq))

data6<-data6 %>%group_by(Age) %>%mutate(sumfreq=sum(Freq))
tmp<-data.frame()
for (i in unique(data6$mutation)){
  for (j in unique(data6$Age)){
  tmp2<-subset(data6, mutation ==i)
  tmp2<-prop.test(x=unique(tmp2$meanfreq),n= unique(tmp2$sumfreq))
  tmp<-rbind(tmp, data.frame(p = tmp2$p.value, mutation = i, Age = j, ymin = tmp2$conf.int[1], ymax = tmp2$conf.int[2]))
  }
}

tmp$p.adj<-p.adjust(tmp$p)
tmp$p.adj<-formatC(tmp$p.adj, digits = 2, format = "e")

#combine mutations of same class





data6$Age<-factor(data6$Age, levels = c("Young","Old"))


dodge<-position_dodge(width=0.7)
data7<-join(data6, tmp)
fig3a<-ggplot(data7, aes(x = mutation, y = meanprop, group = Age,fill = Age, label = p.adj))+
  geom_bar(position = "dodge",col = "black", size = 0.8,width=0.7,stat = "identity",aes(fill  = Age, col = Age))+
  theme_classic()+
  scale_fill_manual(values=c("royalblue","darkorchid4"))+
  ylab("Proportion of each mutation class within cell type")+
  geom_text(y =0.35, check_overlap = T)+ylim(0,0.35)+xlab("Mutation class")+geom_errorbar(position = dodge, width=0.5,aes(ymin=meanprop-seprop, ymax=meanprop+seprop, group = Age))
fig3a
ggsave("220725_Figure_3.pdf", width = 9.5, height = 6.5)








#now do the rest of figure 3####



DDR<- (scan("~/Dropbox (Dropbox @RU)/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/DDR2.txt", what="", sep="\n" ))
testis.list<-ScaleData(testis.list, assay  ="RNA", split.by = "rep")
testis.list@assays[["RNA"]]@scale.data[1:5,1:5]


median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}


tmp<-data.frame()
for ( i in c("Old", "Young")){
  for (j in unique(testis.list$ident)) {
    tmp2<-subset(testis.list, Age == i & ident == j )@assays$RNA@scale.data
    tmp3<-data.frame(gene = rownames(tmp2), exp = rowMeans(tmp2), Age = i, celltype = j)
    tmp<-rbind(tmp, tmp3)
  }
}

tmp$celltype<-factor(x = tmp$celltype, levels=c("Hub cells", "Cyst cells", "Epithelial cells","Accessory gland", "GSC, Early spermatogonia", "Late spermatogonia", "Early spermatocytes", "Late spermatocytes", "Early spermatids", "Late spermatids"))
tmp<-subset(tmp, !celltype %in% c("Hub cells", "Epithelial cells", "Cyst cells","Accessory gland"))
tmp$DDR<-ifelse(tmp$gene %in%  DDR, "Genome\nmaintenance\ngenes", "Other\ngenes")
ggplot(tmp, aes(x= celltype, y = exp, fill = DDR)) + ylim(-1.5,1.5)+
  theme_classic()+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Age, nrow = 2)+theme(axis.text=element_text(size=10))



#add p values
library(rstatix)
stat.test <- tmp %>%
  group_by(celltype,DDR) %>%
  wilcox_test(exp ~ Age) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>% add_y_position(fun ="max")
stat.test$y.position[]<-1
stat.test$p.adj<-formatC(stat.test$p.adj, digits = 2, format = "e")


tmp$Age<-factor(tmp$Age, levels=c("Young","Old"))
#best one
fig3b<-ggplot(tmp, aes(x= DDR, y = exp)) + ylim(-1.5,1.5)+
  theme_classic()+geom_violin(scale = "width",aes(fill = Age))+labs(fill = "Fly age", y = "Scaled expression",x= "Gene type")+
  facet_wrap(~celltype, nrow = 1)+theme(axis.text=element_text(size=10))+stat_summary(aes( x=DDR, y=exp, fill = Age),
                                                                                      fun.data="median_IQR", position=position_dodge(width=0.9), col="white", size=0.5 )+geom_text(check_overlap = T,aes(x, y, label=lab),data=data.frame(x = rep(c("Genome\nmaintenance\ngenes","Other\ngenes"),6), y=Inf, celltype = stat.test$celltype, lab=stat.test$p.adj, yy=letters[1:3]), vjust=1)+scale_fill_manual(values=c("royalblue","darkorchid4"))+xlab("")
plot_grid(fig3a, fig3b, labels = c("A","B"),nrow = 2)
#ggsave("220710_Figure_3.pdf", width = 11.5, height = 9)


Idents(testis.list)<-"ident"
#supplemental figure 1:
p1<-DotPlot(subset(testis.list, Age == "Young"), features = c("p-cup","aub","Prosalpha6T", "His2Av", "ProtB", "Dpy-30L2", "fzo","twe","dlg1","MtnA", "Sems", "Fas3"))+theme(axis.text.x=element_text(angle = 90))
p2<-DotPlot(subset(testis.list, Age == "Old"), features = c("p-cup","aub","Prosalpha6T", "His2Av", "ProtB", "Dpy-30L2", "fzo","twe","dlg1","MtnA", "Sems", "Fas3"))+theme(axis.text.x=element_text(angle = 90))
plot_grid(p1,p2, labels =c("A","B"), nrow = 2)
ggsave("Supplemental_Figure_1.pdf", width=13, height=8)



#supplemental figure 4####
#SNPS vs expression for every replicate

#read in gtf
#bedtools with exons
library(stringr)
gtf<-read.delim("~/Downloads/dmel-all-r6.34.gtf.gz", header = F)
gtf$gene<-str_split(gtf$V9, pattern = "gene_symbol ",simplify = T )[,2]
gtf$gene<-str_split(gtf$gene, pattern = ";",simplify = T )[,1]
gtf<-subset(gtf, V3 =="gene")

#get overlaps
library("GenomicRanges")
library("GRanges")

data6<-data.frame()
tmp3<-unique(data5[,c("chr","pos","ref","alt", "rep")])
for (i in 1:nrow(tmp3)){
  tmp<-unique(tmp3[i,])
  tmp2<-subset(gtf, V1==tmp$chr & V4 < tmp$pos & V5 > tmp$pos )
  if (nrow(tmp2)>0){
    data6<-rbind(data6, data.frame(chr = tmp$chr, pos =tmp$pos, ref = tmp$ref, alt = tmp$alt, gene = tmp2$Gene.name, rep = tmp$rep ))
  }
}





Avgexp<-AverageExpression(testis.list )
Avgexp<-data.frame(Avgexp$RNA)
Avgexp$gene<-rownames(Avgexp)

data7<-join( data6, Avgexp, by = "gene")
data7<-na.omit(data7)
data7$Avgexp<-rowMeans(data7[,7:16])
data7<-data7 %>% group_by(gene) %>% mutate(n=n())
#group by Age, bin by expression, check numbers


data7$rep<-gsub(data7$rep,pattern = "Old3a",replacement = "Old2")
data7$rep<-gsub(data7$rep,pattern = "Old3b",replacement = "Old3")
data7$rep<-gsub(data7$rep,pattern = "Rep1",replacement = "Young1")
data7$rep<-gsub(data7$rep,pattern = "Rep2",replacement = "Young2")
data7$rep<-gsub(data7$rep,pattern = "Younglib3",replacement = "Young3")

data7 %>%ggplot(aes(x=Avgexp, y = n))+geom_point()+facet_wrap(~rep)+theme_classic()+geom_smooth()+ylim(0,50)

data7$bin<-ifelse(data7$Avgexp>mean(data7$Avgexp), "High expression","Low expression")
data7 %>%ggplot(aes(x=bin, y = n, fill = bin))+geom_violin()+geom_boxplot(width=0.1, fill = "white")+facet_wrap(~rep)+theme_classic()+geom_smooth()+ylim(0,50)

#check wilcox
tmp2<-data.frame()
for (i in unique(data7$rep)){
  tmp<-subset(data7, rep ==i)
  p<-wilcox.test(subset(tmp, bin =="High expression")$Avgexp, subset(tmp, bin =="Low expression")$Avgexp)$p.value
  tmp2<-rbind(tmp2, data.frame(p = p, rep = i))
}
tmp2$p<-formatC(tmp2$p, digits = 2)
data7<-join(tmp2, data7, by = "rep")
data7 %>%ggplot(aes(x=bin, y = n, fill = bin))+geom_violin()+
  geom_boxplot(outlier.shape =NA,width=0.1, fill = "white")+
  geom_text(y = 25,check_overlap = T, stat = "identity", aes(x=1.5,label = p))+
  facet_wrap(~rep)+
  theme_classic()+theme(legend.position = 'none')+scale_fill_manual(values=c("royalblue","darkorchid4"))
ggsave("Supplemental_Figure_3.pdf")

#### Supp 3####
