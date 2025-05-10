# BiocManager::install("Seurat")
# 导入数据
pbmc.data <- Read10X(data.dir = "D:/生信/PBMC3Kdata")
# 初始化seurat数据
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)
#查看 pbmc
pbmc  
# 将QC指标可视化，并使用这些指标过滤细胞
#使用PercentageFeatureSet函数计算线粒体QC指标
pbmc[['percent.mt']]<-PercentageFeatureSet(pbmc,pattern = "^MT-")
# 使用violin plot可视化  QC指标，并使用这些指标过滤细胞
VlnPlot(pbmc,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
# FeatureScatter 通常用于可视化两个特征之间的关系
plot1<-FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1+plot2
# 过滤基因数大于2500或少于200的细胞，过滤线粒体基因表达占总表达量>5%的细胞
pbmc<-subset(pbmc,subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)
# 查看过滤后的seurat数据
pbmc 
# 数据规范化
pbmc<-NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
# 识别高度可变的特征（特征选择）
pbmc<-FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)
# 确定高表达的前十个基因
top10<-head(VariableFeatures(pbmc),10)
plot1<-VariableFeaturePlot(pbmc)
plot2<-LabelPoints(plot = plot1,points = top10,repel = TRUE)
plot1+plot2
# 缩放数据
all.genes<-rownames(pbmc)
pbmc<-ScaleData(pbmc,features = all.genes)
# 线性降维
pbmc<-RunPCA(pbmc,features = VariableFeatures(object=pbmc))
print(pbmc[["pca"]],dims = 1:5,nfeatures = 5)
# Seurat提供可视化细胞和定义PCA,包括功能的几种有用的方法ViziDimReduction(),DimPlot()和DimHeatmap()
VizDimLoadings(pbmc,dims = 1:2,reduction = "pca")
DimPlot(pbmc,reduction = "pca")
DimHeatmap(pbmc,dims = 1,cells = 500,balanced = TRUE)
DimHeatmap(pbmc,dims = 1:15,cells=500,balanced = TRUE)
# 确定数据集的维度
pbmc<-JackStraw(pbmc,num.replicate = 100)
pbmc<-ScoreJackStraw(pbmc,dims = 1:20)
# 可视化处理
JackStrawPlot(pbmc,dims = 1:15)
ElbowPlot(pbmc)
# 聚类细胞
# 建立SNN图，并基于其局部领域中的共享重叠细化任意两个单元之间的边权重
pbmc<-FindNeighbors(pbmc,dims = 1:10)
#对细胞进行聚类，应用模块化优化技术，Louvain算法（默认）或SLM
# 以迭代方式将细胞分组在一起，目标是优化标准模块化函数
pbmc<-FindClusters(pbmc,resolution = 0.5)
head(Idents(pbmc),5)
# 运行非线性降维（UMAP/tSNE）
pbmc<-RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = "umap")
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 寻找差异表达的特征（簇生物标志物）
# findmarkers为所有集群自动执行此过程，也可以测试集群组之间的对比，或针对所有单元格进行测试
# 默认情况下，ident.1与所有其他细胞相比，他识别单个簇的阳性和阴性标记。
# min.pct参数要求在两组细胞中的任何一组中以最小百分比检测到一个特征
# 而 thresh.test 参数要求一个特征在两组之间差异表达（平均）一定量
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
# 寻找cluster2的所有markers
cluster2.markers<-FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25)
head(cluster2.markers,n=5)

# 寻找cluster5中与cluster0和cluster3n不同的所有markers
cluster5.markers<-FindMarkers(pbmc,ident.1 = 5,ident.2=c(0,3),min.pct = 0.25)
head(cluster5.markers,n=5)

# 找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞 
pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2,order_by = avg_log2FC)
# 还有用于可视化标记表达的工具，Vlnplot(显示跨集群的表达概率分布)和FeaturePlot()(在tSNE或PCA图上可视化特征表达)是最常用的可视化，建议探索RidgePlot(),CellScatter（），和DotPlot()作为查看数据集的其他方法
VlnPlot(pbmc,features = c("MS4A1","CD79A"))
# 可以绘制行数据
VlnPlot(pbmc,features = c("NKG7","PF4"),slot = "counts",log = TRUE)
# 多个特征图
FeaturePlot(pbmc,features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
#DoHeatmap()为给定的细胞和特征生成一个表达热图。在这种情况下，我们为每个集群绘制前 20 个标记（或所有标记，如果小于 20）。
pbmc.markers%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC)->top10
DoHeatmap(pbmc,features = top10$gene)+NoLegend()
# 将细胞类型标识分配给集群
new.cluster.ids<-c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                   "NK", "DC", "Platelet")
names(new.cluster.ids)<-levels(pbmc)
pbmc<-RenameIdents(pbmc,new.cluster.ids)
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
saveRDS(pbmc,file = "./pbmc3k_final.rds")
