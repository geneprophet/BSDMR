# mouse_met <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/mouse_8",min_coverage=4)
# mouse_anno <- read_annotation("/data40T/kanghe/R/Mouse.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/Human.gtf")
# human_met <- read_methylation_report("/data40T/kanghe/WGBS/human/human_10_CX_report.txt",min_coverage=4)
# human_region <- cluster_sites_to_region(methylation = human_met,annotation = human_anno,is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_obj <- create_region_object(human_met, human_region)
# human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_obj <- create_region_object(mouse_met, mouse_region)
# mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100,no_cores=50)
# save(mouse_fit_mean,file = '/pnas/zhangz_group/lizhao/mouse_mean.Rdata')
# similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,subjectProfiles=mouse_fit_profiles)

# mouse_region2 <- cluster_sites_to_region(methylation = mouse_met,annotation = mouse_anno,is_parallel = TRUE)
# filePath <- system.file("extdata", "mm10ToHg38.over.chain", package = "BSDMR")
# human_region2 <- change_genomic_coordinate(mouse_region2,filePath,human_anno)
# human_obj2 <- create_region_object(human_met, human_region2)
# mouse_obj2 <- create_region_object(mouse_met, mouse_region2)
# mouse_fit_profiles2 <- infer_profiles_vb(X = mouse_obj2$met,basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100,no_cores = 20)
# mouse_fit_mean2 <- infer_profiles_vb(X = mouse_obj2$met,basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_profiles2 <- infer_profiles_vb(X = human_obj2$met,basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean2 <- infer_profiles_vb(X = human_obj2$met, basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
###########brain
# mouse_anno <- read_annotation("/data40T/kanghe/R/Mouse.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/Human.gtf")
# human_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/human/brain/human_12_brain",min_coverage=4)
# mouse_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/brain/mouse_6_brain",min_coverage=4)
# human_region_brain <- cluster_sites_to_region(methylation = human_met_brain, annotation = human_anno, is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region_brain <- change_genomic_coordinate(human_region_brain,filePath,mouse_anno)
# human_obj_brain <- create_region_object(human_met_brain, human_region_brain)
# mouse_obj_brain <- create_region_object(mouse_met_brain, mouse_region_brain)
#

# which(similarity<0.1)
# i = 35165                                                                                                                                                                                
# plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#                     obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
# plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#                     obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))
# 

# ############可视化
# a=width(human_region[which(similarity<0.2)])
# species = rep("human_DMR",NROW(a))
# a = data.frame(width=a,species=species)
# b=width(mouse_region[which(similarity<0.2)])
# species = rep("mouse_DMR",NROW(b))
# b = data.frame(width=b,species = species)
# c=width(human_region[which(similarity>0.9)])
# species = rep("human_IMR",NROW(c))
# c = data.frame(width=c,species=species)
# d=width(mouse_region[which(similarity>0.9)])
# species = rep("mouse_IMR",NROW(d))
# d = data.frame(width=d,species = species)
# e = rbind(a,b,c,d)
# 
# ggplot(e, aes(x = factor(species), y = width, fill = factor(species))) +
#   geom_boxplot(notch = TRUE) +
#   scale_fill_brewer(palette = "Pastel2") +
#   xlab(label = "Type") +
#   ylab(label = "Region Width") +
#   labs(fill = "Species")
# 
# #
# a = unique(human_region[which(similarity<0.2)]$id)
# b = unique(human_region[which(similarity>0.9)]$id)
# veen <- venn.diagram(list(DMR=a,IMR=b),
#                      filename = 'DMR_IMR.tiff',
#                      height = 3000,
#                      width = 3000,
#                      resolution = 500,
#                      col = "transparent",
#                      fill = c("red","blue"),
#                      alpha = 0.5,
#                      label.col = c("darkred","white","darkblue"),
#                      cex = 2.5,
#                      fontfamily = "serif",
#                      fontface = "bold",
#                      cat.default.pos = "outer",#设置标签在圆外面
#                      cat.cex = 2,#外标签的字体大小
#                      cat.fontfamily = "serif",
#                      cat.dist = c(-0.05,-0.02)
#                      )



# a=Homo_sapiens_TF$Ensembl
# b = unique(human_region[which(similarity<0.2)]$id)
# c = unique(human_region[which(similarity>0.9)]$id)
# veen <- venn.diagram(list(TF=a,DMR=b,IMR=c),
#                      filename = 'TF_IMR_DMR.tiff',
#                      height = 3000,
#                      width = 3000,
#                      resolution = 500,
#                      col = "transparent",
#                      fill = c("red","blue","green"),
#                      alpha = 0.5,
#                      #label.col = c("darkred","white","darkblue"),
#                      cex = 2.5,
#                      fontfamily = "serif",
#                      fontface = "bold",
#                      cat.default.pos = "outer",#设置标签在圆外面
#                      cat.cex = 2,#外标签的字体大小
#                      cat.fontfamily = "serif",
#                      #cat.dist = c(-0.05,-0.02)
# )

#

# ###GO_KEGG analysis
# #
# # library(org.Hs.eg.db)
# a = unique(human_region[which(similarity<0.2)]$id)
# b = unique(human_region[which(similarity>0.9)]$id)
# c <- intersect(a,b)
# #gene <- setdiff(a,c)
# gene <- setdiff(b,c)
# #gene <- c
# allGene <- unique(human_region$id)
# geneList <- bitr(allGene,fromType = "ENSEMBL",
#                  toType = "ENTREZID",
#                  OrgDb = org.Hs.eg.db)
# gene_df <- bitr(gene, fromType = "ENSEMBL",
#                 toType = c("ENTREZID", "SYMBOL"),
#                 OrgDb = org.Hs.eg.db)
# head(gene_df)
# ego <- enrichGO(gene          = gene_df$ENTREZID,
#                 universe      = geneList$ENTREZID,
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05)
# head(summary(ego))
# dotplot(ego, showCategory=30)
# 
# kk <- enrichKEGG(gene          = gene_df$ENTREZID,
#                  universe      = geneList$ENTREZID,
#                organism = 'hsa',
#                pvalueCutoff=0.05,
#                pAdjustMethod="BH",
#                qvalueCutoff=0.1)
# head(kk)
# dotplot(kk, showCategory=30)

#3k liver human and mouse
# mouse_anno <- read_annotation("/data40T/kanghe/R/MOUSE_3K.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/HUMAN_3K.gtf")
# mouse_met <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/mouse_8",min_coverage=4)
# human_met <- read_methylation_report("/data40T/kanghe/WGBS/human/human_10_CX_report.txt",min_coverage=4)
# human_region <- cluster_sites_to_region(methylation = human_met,annotation = human_anno,is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_obj <- create_region_object(human_met, human_region)
# human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_obj <- create_region_object(mouse_met, mouse_region)
# mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,subjectProfiles=mouse_fit_profiles)
# 
# which(similarity<0.2)
# i = 140
# plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#                     obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
# plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#                     obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))

############################################
#3k brain human and mouse
# mouse_anno <- read_annotation("/data40T/kanghe/R/MOUSE_3K.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/HUMAN_3K.gtf")
# human_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/human/brain/human_12_brain",min_coverage=4)
# mouse_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/brain/mouse_6_brain",min_coverage=4)
# human_region_brain <- cluster_sites_to_region(methylation = human_met_brain, annotation = human_anno, is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region_brain <- change_genomic_coordinate(human_region_brain,filePath,mouse_anno)
# human_obj_brain <- create_region_object(human_met_brain, human_region_brain)
# mouse_obj_brain <- create_region_object(mouse_met_brain, mouse_region_brain)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_fit_profiles_brain <- infer_profiles_vb(X = human_obj_brain$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean_brain <- infer_profiles_vb(X = human_obj_brain$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_fit_profiles_brain <- infer_profiles_vb(X = mouse_obj_brain$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean_brain <- infer_profiles_vb(X = mouse_obj_brain$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# similarity_brain <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles_brain,subjectProfiles=mouse_fit_profiles_brain)



# #ChipSeeker Annotation
load(file = '../liver_3k.Rdata')
h_dmr = human_region[which(similarity<0.2)]
m_dmr = mouse_region[which(similarity<0.2)]
length(which(h_dmr$annotation=="gene"))
length(which(h_dmr$annotation=="promoter"))
length(which(h_dmr$annotation=="downstream"))
h_imr = human_region[which(similarity>0.9)]
m_imr = mouse_region[which(similarity>0.9)]
length(which(h_imr$annotation=="gene"))
length(which(h_imr$annotation=="promoter"))
length(which(h_imr$annotation=="downstream"))
covplot(h_dmr)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
h_dmr_anno = annotatePeak(h_dmr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
plotAnnoPie(h_dmr_anno)
vennpie(h_dmr_anno)
upsetplot(h_dmr_anno)

covplot(h_imr)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
h_imr_anno = annotatePeak(h_imr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
plotAnnoPie(h_imr_anno)
vennpie(h_imr_anno)
upsetplot(h_imr_anno)


GRCh38_ccREs <- read_delim("D:/R/R_project/GRCh38-ccREs.bed",
                            "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
hg38CRE <- GenomicRanges::GRanges(
  seqnames = GRCh38_ccREs$X1,
  ranges = IRanges::IRanges(start = GRCh38_ccREs$X2, end = GRCh38_ccREs$X3),
  type = GRCh38_ccREs$X6
)
#查看CRE在基于组上的分布情况
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
h_CRE_anno <- annotatePeak(hg38CRE,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
plotAnnoPie(h_CRE_anno)


HUMAN_3K <- read_delim("D:/R/R_project/HUMAN_3K.gtf",
                       "\t", escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)
hg38 <- GenomicRanges::GRanges(
  seqnames = paste0("chr",HUMAN_3K$X1),
  ranges = IRanges::IRanges(start = HUMAN_3K$X4, end = HUMAN_3K$X5),
  EnsemblID = HUMAN_3K$X2
)

##########################################################################

#去除重复的DMR
length(unique(paste0(h_dmr$center,width(h_dmr))))
length(h_dmr)
h_dmr = h_dmr[!duplicated(h_dmr$center)]

##DMR所在基因的DMR和CRE的富集分析
overlaps1 <- GenomicRanges::findOverlaps(query = hg38,subject = hg38CRE,ignore.strand = T)
overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_dmr,ignore.strand = T)
length(unique(queryHits(overlaps1)))
length(unique(queryHits(overlaps2)))


kkDMR <- list()
for (i in 1:length(unique(queryHits(overlaps2)))) {
  #当前包含DMR基因
  gene <- hg38[unique(queryHits(overlaps2))[i]]
  #当前基因的所有DMR
  DMR <- h_dmr[subjectHits(overlaps2[which(queryHits(overlaps2)==unique(queryHits(overlaps2))[i])])]
  DMR <- GenomicRanges::intersect(DMR,gene,ignore.strand=TRUE)
  
  #当前基因的所有CRE
  CRE <- hg38CRE[subjectHits(overlaps1[which(queryHits(overlaps1)==unique(queryHits(overlaps2))[i])])]
  CRE <- GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE)
  
  a = ifelse(length(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))),0)
  
  b = sum(width(DMR)) - a 
  
  c = sum(width(CRE)) - a 
  
  d = width(gene) - (a+b+c)
  
  data <- data.frame(YesDMR=c(a,b),NoDMR=c(c,d))
  rownames(data)=c("YesCRE","NoCRE")
  fisher <- fisher.test(data)
  kkDMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,DMR=DMR,CRE=CRE)

}

#查看p_val和OR值
DMR_p_val <- vector(mode = "numeric",length = 0)
DMR_odds_ratio <- vector(mode = "numeric",length = 0)

for (i in 1:length(kkDMR)) {
  DMR_p_val <- c(DMR_p_val,kkDMR[[i]]$p_val)
  DMR_odds_ratio <- c(DMR_odds_ratio,kkDMR[[i]]$odds_ratio)
}

length(DMR_p_val[which(DMR_p_val<0.01)])
DMR_odds_ratio <- DMR_odds_ratio[which(DMR_p_val<0.01)]
length(DMR_odds_ratio[which(DMR_odds_ratio<1)])
length(DMR_odds_ratio[which(DMR_odds_ratio>1)])
length(kkDMR)


##########################################################################
#去除重复的IMR
length(unique(paste0(h_imr$center,width(h_imr))))
length(h_imr)
h_imr = h_imr[!duplicated(h_imr$center)]

##IMR所在基因的IMR和CRE的富集分析
overlaps1 <- GenomicRanges::findOverlaps(query = hg38,subject = hg38CRE,ignore.strand = T)
overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_imr,ignore.strand = T)
length(unique(queryHits(overlaps1)))
length(unique(queryHits(overlaps2)))

kkIMR <- list()
for (i in 1:length(unique(queryHits(overlaps2)))) {
  #当前包含IMR基因
  #gene <- width(hg38[unique(queryHits(overlaps2))[i]])
  gene <- hg38[unique(queryHits(overlaps2))[i]]
  #当前基因的所有IMR
  IMR <- h_imr[subjectHits(overlaps2[which(queryHits(overlaps2)==unique(queryHits(overlaps2))[i])])]
  IMR <- GenomicRanges::intersect(IMR,gene,ignore.strand=TRUE)
  
  #当前基因的所有CRE
  CRE <- hg38CRE[subjectHits(overlaps1[which(queryHits(overlaps1)==unique(queryHits(overlaps2))[i])])]
  CRE <- GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE)
  
  
  a = ifelse(length(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))),0)
  
  b = sum(width(GenomicRanges::intersect(IMR,gene,ignore.strand=TRUE))) - a 
  
  c = sum(width(GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE))) - a 
  
  d = width(gene) - (a+b+c)

  data <- data.frame(YesIMR=c(a,b),NoIMR=c(c,d))
  rownames(data)=c("YesCRE","NoCRE")
  fisher <- fisher.test(data)
  kkIMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,IMR=IMR,CRE=CRE)

}

#查看p_val和OR值
IMR_p_val <- vector(mode = "numeric",length = 0)
IMR_odds_ratio <- vector(mode = "numeric",length = 0)

for (i in 1:length(kkIMR)) {
  IMR_p_val <- c(IMR_p_val,kkIMR[[i]]$p_val)
  IMR_odds_ratio <- c(IMR_odds_ratio,kkIMR[[i]]$odds_ratio)
}

length(IMR_p_val[which(IMR_p_val<0.01)])
IMR_odds_ratio <- IMR_odds_ratio[which(IMR_p_val<0.01)]
length(IMR_odds_ratio[which(IMR_odds_ratio<1)])
length(IMR_odds_ratio[which(IMR_odds_ratio>1)])
length(kkIMR)

#整体看DMR\IMR和CRE的富集分析
sum(width(GenomicRanges::intersect(h_dmr,hg38CRE,ignore.strand=T)))
sum(width(h_dmr))
sum(width(hg38CRE))
sum(width(hg38))
fisher.test(matrix(c(1287416,3683312,252960490,797919039),2,2))

sum(width(GenomicRanges::intersect(h_imr,hg38CRE,ignore.strand=T)))
sum(width(h_imr))
sum(width(hg38CRE))
sum(width(hg38))
fisher.test(matrix(c(3493537,8694956,250754369,792907395),2,2))



#overlap of DMR and IMR , overlaps gene of DMR_OR<1 and IMR_OR<1
overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_imr,ignore.strand = T)
EnsemblID_IMR <- unique(hg38[queryHits(overlaps2)]$EnsemblID)
overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_dmr,ignore.strand = T)
EnsemblID_DMR <- unique(hg38[queryHits(overlaps2)]$EnsemblID)
length(EnsemblID_DMR)
length(EnsemblID_IMR)
length(base::intersect(EnsemblID_DMR,EnsemblID_IMR))



OR_less1 <- EnsemblID_DMR[which(odds_ratio<1)]
#查看交集
length(intersect(EnsemblID_IMR,OR_less1))
#查看差异
length(setdiff(EnsemblID_IMR,OR_less1))
length(OR_less1)
#TODO:// 查看OR大于1的DMR所在基因的功能富集分析

#查看OR值大于1和小于1的DMR的长度分布，有没有可能是长度过长导致的CRE富集
OR_more1 <- EnsemblID_DMR[which(odds_ratio>1)]
length(OR_more1)



#####################################################################################
#查看脑组织的情况
load('../brain_3k.Rdata')
