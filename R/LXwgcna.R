LXwgcna <- function(data_file,expriment_ifo,key_phenotype){

#-----------安装和导入相关R包（用sapply函数）
pack_install <- function(){

  installed_packages <- data.frame(installed.packages()) #查看已安装的R包

  comman_pack <- c("BiocManager","purrr","openxlsx","dplyr","magrittr","ggplot2","tidyr","ggplotify") #普通R包 https://cran.rstudio.com/
  Bio_pack <- c("WGCNA","impute","preprocessCore","limma","Biostrings","Biobase","GO.db") # 生物医学相关R包 https://bioconductor.org/packages/
  #--"Biostrings","Biobase"是"GO.db"依赖包

  # 判断是否已安装，如果没装，则提取出来
  comman_install <- comman_pack[!comman_pack %in% installed_packages$Package]
  Bio_install <- Bio_pack[!Bio_pack %in% installed_packages$Package]

  # 用sapply函数批量安装R包
  comman_fun <- function(i){install.packages(i,update=F,ask=F) }
  sapply(comman_install,comman_fun,simplify = T)

  Bio_fun <- function(i){BiocManager::install(i,update=F,ask=F)}
  sapply(Bio_install,Bio_fun,simplify = T)

  # 批量library
  library_fun <- function(i){library(i,character.only = T)}
  sapply(c(comman_pack,Bio_pack),library_fun,simplify = T)

}

pack_install()

#------处理表达数据------------------------------------------------------------
if (!dir.exists("analysis results"))
    dir.create("analysis results")

exp_df0 <- read.xlsx(data_file)
exp_df <- limma::avereps (exp_df0[,-1],ID = exp_df0[,1]) %>% #对代谢代ID去重，同时取平均值
  data.frame() %>% na.omit()

datExpr <- dplyr::filter(exp_df,rowSums(exp_df)>0) %>% #去掉表达值都是0的行
  t() %>% data.frame() #行和列转置

# 聚类所有样本，观察是否有离群值或异常值。
while (!is.null(dev.list()))  dev.off()
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/1.Sample clustering.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

while (!is.null(dev.list()))   dev.off()
# -----------------载入表型数据-----------------
traitData = read.xlsx(expriment_ifo,rowNames = T) %>%
  t() %>% data.frame()

#可视化表型数据与基因表达量数据的联系，重构样本聚类树
# 颜色越深，代表这个表型数据与这个样本的基因表达量关系越密切
traitColors = numbers2colors(traitData, signed = FALSE) #用颜色代表关联度
png(filename = "analysis results/2.Sample dendrogram and trait heatmap.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")
while (!is.null(dev.list()))  dev.off()

# 构建表达网络
# 构建表达网络是WGCNA分析中最为关键的一步，是否构建成功、
#是否构建正确对后期模块的划分和关联表型数据筛选核心基因至关重要。
#挑选软阈值是构建网络拓扑分析的关键，选择软阈值是基于近无尺度拓扑标准的。

###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
png(filename = "analysis results/3.scale independencep.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red") #可以修改

while (!is.null(dev.list()))  dev.off()
###平均连通性与power值散点图
png(filename = "analysis results/4.Mean Connectivity.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

while (!is.null(dev.list()))  dev.off()

###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
adjacency = adjacency(datExpr, power = softPower)
softPower

#无尺度网络的拓扑关系
#pdf(file="3_softConnectivity.pdf", width=9, height=5)
k <- softConnectivity(datE=datExpr,power=softPower)
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
#dev.off()

###TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

###基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Metabolite clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
#dev.off()

###动态剪切模块识别
minModuleSize = 10      #模块数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
png(filename = "analysis results/5.Metabolite dendrogram and module colors.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plotDendroAndColors(geneTree, dynamicColors, "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolite dendrogram and module colors")

while (!is.null(dev.list()))   dev.off()

###对模块进行聚类,找出相似模块聚类
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
#pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25  #剪切高度可修改
abline(h=MEDissThres, col = "red")
#dev.off()

###相似模块合并
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
png(filename = "analysis results/6.Merged metabolite dendrogram and module colors.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plotDendroAndColors(geneTree, mergedColors,"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolite dendrogram and module colors")

while (!is.null(dev.list()))  dev.off()

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

###模块与性状数据热图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, traitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#pdf(file="8_Module_trait.pdf", width=6, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

png(filename = "analysis results/7.Module-trait relationships.png",
    width=1000, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

while (!is.null(dev.list()))  dev.off()

###计算MM和GS值
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(traitData)
geneTraitSignificance = as.data.frame(cor(datExpr, traitData, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###输出模块重要性的图形
y=traitData[,1]
GS1=as.numeric(cor(y, datExpr, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, mergedColors, mean, na.rm=T)
#pdf(file="9_GeneSignificance.pdf", width=11, height=7)
sizeGrWindow(10, 8)
png(filename = "analysis results/8.Metabolite significance across modules.png",
    width=1000, height=600,units = "px",res = 100)
par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)
plotModuleSignificance(geneSignificance = GeneSignificance, colors = mergedColors,boxplot=TRUE,
                       main = "Metabolite significance across modules,",ylab = "Metabolite Significance")
while (!is.null(dev.list()))  dev.off()

###批量输出性状和模块散点图
traitColumn=match(key_phenotype,traitNames)

for (module in modNames){
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  if (nrow(geneModuleMembership[moduleGenes,]) > 1){
    outPng=paste("analysis results/9_", key_phenotype, "_", module,".png",sep="")
   # pdf(file=outPdf,width=7,height=7)
   # par(mfrow = c(1,1))
    png(filename = outPng, width=1000, height=600,units = "px",res = 100)
    par(mai = c(1, 2, 1, 1)) # 边距 c(bottom, left, top, right)

    verboseScatterplot(abs(as.numeric(geneModuleMembership[moduleGenes, column])),
                       abs(as.numeric(geneTraitSignificance[moduleGenes, traitColumn])),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Metabolite significance for ",key_phenotype),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
                       abline(v=0.8,h=0.5,col="red")
    while (!is.null(dev.list()))  dev.off()

    }
}


###输出GS_MM数据
probes = colnames(datExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.xlsx(geneInfo, "analysis results/GS_MM.xlsx")


###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, paste0("analysis results/","module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###输出每个模块的核心基因
geneSigFilter=0.5         #基因重要性的过滤条件
moduleSigFilter=0.8       #基因与模块相关性的过滤条件
datMM=cbind(geneModuleMembership, geneTraitSignificance)
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
for(mmi in colnames(datMM)[1:(ncol(datMM)-2)]){
  dataMM2=datMM[abs(datMM[,mmi])>moduleSigFilter,]
  write.table(row.names(dataMM2), file =paste0("analysis results/","hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

print("----------------------------------------------------------------------")
print("The results of WGCNA can be found in the folder of <analysis results> ")
print("----------------------------------------------------------------------")
}
