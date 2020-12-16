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
# which(similarity<0.1)
# i = 297641               
# plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#                     obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
# plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#                     obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))

# if(F){
#   load('../liver_3k.Rdata')
#   ############可视化
#   a=width(human_region[which(similarity<0.2)])
#   species = rep("human_DMR",NROW(a))
#   a = data.frame(width=a,species=species)
#   b=width(mouse_region[which(similarity<0.2)])
#   species = rep("mouse_DMR",NROW(b))
#   b = data.frame(width=b,species = species)
#   c=width(human_region[which(similarity>0.9)])
#   species = rep("human_IMR",NROW(c))
#   c = data.frame(width=c,species=species)
#   d=width(mouse_region[which(similarity>0.9)])
#   species = rep("mouse_IMR",NROW(d))
#   d = data.frame(width=d,species = species)
#   e = rbind(a,b,c,d)
#   
#   ggplot(e, aes(x = factor(species), y = width, fill = factor(species))) +
#     geom_boxplot(notch = TRUE) +
#     scale_fill_brewer(palette = "Pastel2") +
#     xlab(label = "Type") +
#     ylab(label = "Region Width") +
#     labs(fill = "Species")
#   
#   #
#   a = unique(human_region[which(similarity<0.2)]$id)
#   b = unique(human_region[which(similarity>0.9)]$id)
#   veen <- venn.diagram(list(DMR=a,IMR=b),
#                        filename = 'DMR_IMR.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue"),
#                        alpha = 0.5,
#                        label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        cat.dist = c(-0.04,-0.1)
#   )
#   
#   a=Homo_sapiens_TF$Ensembl
#   b = unique(human_region[which(similarity<0.2)]$id)
#   c = unique(human_region[which(similarity>0.9)]$id)
#   veen <- venn.diagram(list(TF=a,DMR=b,IMR=c),
#                        filename = 'TF_IMR_DMR.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue","green"),
#                        alpha = 0.5,
#                        #label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        #cat.dist = c(-0.05,-0.02)
#   )
#   
#   ###GO_KEGG analysis
#   #
#   library(org.Hs.eg.db)
#   a = unique(human_region[which(similarity<0.2)]$id)
#   b = unique(human_region[which(similarity>0.9)]$id)
#   c <- intersect(a,b)
#   gene <- setdiff(a,c)
#   #gene <- setdiff(b,c)
#   #gene <- c
#   allGene <- unique(human_region$id)
#   geneList <- bitr(allGene,fromType = "ENSEMBL",
#                    toType = "ENTREZID",
#                    OrgDb = org.Hs.eg.db)
#   gene_df <- bitr(gene, fromType = "ENSEMBL",
#                   toType = c("ENTREZID", "SYMBOL"),
#                   OrgDb = org.Hs.eg.db)
#   head(gene_df)
#   ego <- enrichGO(gene          = gene_df$ENTREZID,
#                   universe      = geneList$ENTREZID,
#                   OrgDb         = org.Hs.eg.db,
#                   ont           = "ALL",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.01,
#                   qvalueCutoff  = 0.05)
#   head(summary(ego))
#   dotplot(ego, showCategory=30)
#   
#   kk <- enrichKEGG(gene          = gene_df$ENTREZID,
#                    universe      = geneList$ENTREZID,
#                    organism = 'hsa',
#                    pvalueCutoff=0.05,
#                    pAdjustMethod="BH",
#                    qvalueCutoff=0.1)
#   head(kk)
#   dotplot(kk, showCategory=30)
#   
#   ############查看GTEx对应组织的top100表达基因
#   library(readr)
#   liverGTEx_Portal <- read_delim("D:/downloads/liverGTEx Portal.txt", 
#                                  "\t", escape_double = FALSE, trim_ws = TRUE)
#   #substr(liverGTEx_Portal$`Gencode Id`,1,15)
#   veen <- venn.diagram(list(DMR=a,TOP100=substr(liverGTEx_Portal$`Gencode Id`,1,15)),
#                        filename = 'liver_DMR_top100.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue"),
#                        alpha = 0.5,
#                        label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        cat.dist = c(-0.05,-0.02)
#   )
#   
# }
# if(F){
#   # #ChipSeeker Annotation
#   load(file = '../liver_3k.Rdata')
  # h_dmr = human_region[which(similarity<0.2)]
  # m_dmr = mouse_region[which(similarity<0.2)]
  # 
  # h_imr = human_region[which(similarity>0.9)]
  # m_imr = mouse_region[which(similarity>0.9)]
  # 
  # covplot(h_dmr)
  # txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  # h_dmr_anno = annotatePeak(h_dmr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
  # plotAnnoPie(h_dmr_anno)
  # vennpie(h_dmr_anno)
  # upsetplot(h_dmr_anno)
  # 
  # covplot(h_imr)
  # txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  # h_imr_anno = annotatePeak(h_imr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
  # plotAnnoPie(h_imr_anno)
  # vennpie(h_imr_anno)
  # upsetplot(h_imr_anno)
#   
# }
# 
# library(readr)
# GRCh38_ccREs <- read_delim("D:/R/R_project/GRCh38-ccREs.bed",
#                             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# hg38CRE <- GenomicRanges::GRanges(
#   seqnames = GRCh38_ccREs$X1,
#   ranges = IRanges::IRanges(start = GRCh38_ccREs$X2, end = GRCh38_ccREs$X3),
#   type = GRCh38_ccREs$X6
# )
# #查看CRE在基于组上的分布情况
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# h_CRE_anno <- annotatePeak(hg38CRE,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
# plotAnnoPie(h_CRE_anno)
# 
# HUMAN_3K <- read_delim("D:/R/R_project/HUMAN_3K.gtf",
#                        "\t", escape_double = FALSE, col_names = FALSE,
#                        trim_ws = TRUE)
# hg38 <- GenomicRanges::GRanges(
#   seqnames = paste0("chr",HUMAN_3K$X1),
#   ranges = IRanges::IRanges(start = HUMAN_3K$X4, end = HUMAN_3K$X5),
#   EnsemblID = HUMAN_3K$X2
# )
# 
# ##########################################################################
# 
# #去除重复的DMR
# length(unique(paste0(h_dmr$center,width(h_dmr))))
# length(h_dmr)
# h_dmr = h_dmr[!duplicated(h_dmr$center)]
# 
# ###################################################
# if(F){
#   DMR_ID <- unique(h_dmr$id)
#   length(DMR_ID)
#   kkDMR <- list()
#   for (i in 1:length(DMR_ID)) {
#     #当前包含DMR基因
#     gene <- hg38[which(hg38$EnsemblID==DMR_ID[i])]
#     #当前基因的所有DMR
#     DMR <- GenomicRanges::intersect(h_dmr,gene,ignore.strand=TRUE)
#     #当前基因的所有CRE
#     CRE <- GenomicRanges::intersect(hg38CRE,gene,ignore.strand=TRUE)
#     a = ifelse(length(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))),0)
#     
#     b = sum(width(DMR)) - a 
#     
#     c = sum(width(CRE)) - a 
#     
#     d = width(gene) - (a+b+c)
#     
#     data <- data.frame(YesDMR=c(a,b),NoDMR=c(c,d))
#     rownames(data)=c("YesCRE","NoCRE")
#     fisher <- fisher.test(data)
#     kkDMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,DMR=DMR,CRE=CRE)
#     
#   }
#   
#   #查看p_val和OR值
#   DMR_p_val <- vector(mode = "numeric",length = 0)
#   DMR_odds_ratio <- vector(mode = "numeric",length = 0)
#   
#   for (i in 1:length(kkDMR)) {
#     DMR_p_val <- c(DMR_p_val,kkDMR[[i]]$p_val)
#     DMR_odds_ratio <- c(DMR_odds_ratio,kkDMR[[i]]$odds_ratio)
#   }
#   
#   length(DMR_p_val[which(DMR_p_val<0.01)])
#   DMR_odds_ratio <- DMR_odds_ratio[which(DMR_p_val<0.01)]
#   length(DMR_odds_ratio[which(DMR_odds_ratio<1)])
#   length(DMR_odds_ratio[which(DMR_odds_ratio>10)])
#   length(kkDMR)
#   DMR_OR_greter_10 <- DMR_ID[which(DMR_odds_ratio>10)]
#   
# }
# 
# if(F){
#   ##DMR所在基因的DMR和CRE的富集分析
#   overlaps1 <- GenomicRanges::findOverlaps(query = hg38,subject = hg38CRE,ignore.strand = T)
#   overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_dmr,ignore.strand = T)
#   length(unique(queryHits(overlaps1)))
#   length(unique(queryHits(overlaps2)))
#   
#   kkDMR <- list()
#   for (i in 1:length(unique(queryHits(overlaps2)))) {
#     #当前包含DMR基因
#     gene <- hg38[unique(queryHits(overlaps2))[i]]
#     #当前基因的所有DMR
#     DMR <- h_dmr[subjectHits(overlaps2[which(queryHits(overlaps2)==unique(queryHits(overlaps2))[i])])]
#     DMR <- GenomicRanges::intersect(DMR,gene,ignore.strand=TRUE)
#     
#     #当前基因的所有CRE
#     CRE <- hg38CRE[subjectHits(overlaps1[which(queryHits(overlaps1)==unique(queryHits(overlaps2))[i])])]
#     CRE <- GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE)
#     
#     a = ifelse(length(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))),0)
#     
#     b = sum(width(DMR)) - a 
#     
#     c = sum(width(CRE)) - a 
#     
#     d = width(gene) - (a+b+c)
#     
#     data <- data.frame(YesDMR=c(a,b),NoDMR=c(c,d))
#     rownames(data)=c("YesCRE","NoCRE")
#     fisher <- fisher.test(data)
#     kkDMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,DMR=DMR,CRE=CRE)
#     
#   }
#   
#   #查看p_val和OR值
#   DMR_p_val <- vector(mode = "numeric",length = 0)
#   DMR_odds_ratio <- vector(mode = "numeric",length = 0)
#   
#   for (i in 1:length(kkDMR)) {
#     DMR_p_val <- c(DMR_p_val,kkDMR[[i]]$p_val)
#     DMR_odds_ratio <- c(DMR_odds_ratio,kkDMR[[i]]$odds_ratio)
#   }
#   
#   length(DMR_p_val[which(DMR_p_val<0.01)])
#   DMR_odds_ratio <- DMR_odds_ratio[which(DMR_p_val<0.01)]
#   length(DMR_odds_ratio[which(DMR_odds_ratio<1)])
#   length(DMR_odds_ratio[which(DMR_odds_ratio>1)])
#   length(kkDMR)
#   
# }
# 
# ##################################################
# 
# if(F){
#   #去除重复的IMR
#   length(unique(paste0(h_imr$center,width(h_imr))))
#   length(h_imr)
#   h_imr = h_imr[!duplicated(h_imr$center)]
#   
#   ##IMR所在基因的IMR和CRE的富集分析
#   overlaps1 <- GenomicRanges::findOverlaps(query = hg38,subject = hg38CRE,ignore.strand = T)
#   overlaps2 <- GenomicRanges::findOverlaps(query = hg38,subject = h_imr,ignore.strand = T)
#   length(unique(queryHits(overlaps1)))
#   length(unique(queryHits(overlaps2)))
#   
#   kkIMR <- list()
#   for (i in 1:length(unique(queryHits(overlaps2)))) {
#     #当前包含IMR基因
#     #gene <- width(hg38[unique(queryHits(overlaps2))[i]])
#     gene <- hg38[unique(queryHits(overlaps2))[i]]
#     #当前基因的所有IMR
#     IMR <- h_imr[subjectHits(overlaps2[which(queryHits(overlaps2)==unique(queryHits(overlaps2))[i])])]
#     IMR <- GenomicRanges::intersect(IMR,gene,ignore.strand=TRUE)
#     
#     #当前基因的所有CRE
#     CRE <- hg38CRE[subjectHits(overlaps1[which(queryHits(overlaps1)==unique(queryHits(overlaps2))[i])])]
#     CRE <- GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE)
#     
#     
#     a = ifelse(length(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))),0)
#     
#     b = sum(width(GenomicRanges::intersect(IMR,gene,ignore.strand=TRUE))) - a 
#     
#     c = sum(width(GenomicRanges::intersect(CRE,gene,ignore.strand=TRUE))) - a 
#     
#     d = width(gene) - (a+b+c)
#     
#     data <- data.frame(YesIMR=c(a,b),NoIMR=c(c,d))
#     rownames(data)=c("YesCRE","NoCRE")
#     fisher <- fisher.test(data)
#     kkIMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,IMR=IMR,CRE=CRE)
#     
#   }
#   
#   #查看p_val和OR值
#   IMR_p_val <- vector(mode = "numeric",length = 0)
#   IMR_odds_ratio <- vector(mode = "numeric",length = 0)
#   
#   for (i in 1:length(kkIMR)) {
#     IMR_p_val <- c(IMR_p_val,kkIMR[[i]]$p_val)
#     IMR_odds_ratio <- c(IMR_odds_ratio,kkIMR[[i]]$odds_ratio)
#   }
#   
#   length(IMR_p_val[which(IMR_p_val<0.01)])
#   IMR_odds_ratio <- IMR_odds_ratio[which(IMR_p_val<0.01)]
#   length(IMR_odds_ratio[which(IMR_odds_ratio<1)])
#   length(IMR_odds_ratio[which(IMR_odds_ratio>1)])
#   length(kkIMR)
# }
# 
# if(F){
#   IMR_ID <- unique(h_imr$id)
#   length(IMR_ID)
#   kkIMR <- list()
#   for (i in 1:length(IMR_ID)) {
#     #当前包含IMR基因
#     gene <- hg38[which(hg38$EnsemblID==IMR_ID[i])]
#     #当前基因的所有IMR
#     IMR <- GenomicRanges::intersect(h_imr,gene,ignore.strand=TRUE)
#     #当前基因的所有CRE
#     CRE <- GenomicRanges::intersect(hg38CRE,gene,ignore.strand=TRUE)
#     a = ifelse(length(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(IMR,CRE,ignore.strand=TRUE))),0)
#     
#     b = sum(width(IMR)) - a 
#     
#     c = sum(width(CRE)) - a 
#     
#     d = width(gene) - (a+b+c)
#     
#     data <- data.frame(YesIMR=c(a,b),NoIMR=c(c,d))
#     rownames(data)=c("YesCRE","NoCRE")
#     fisher <- fisher.test(data)
#     kkIMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,IMR=IMR,CRE=CRE)
#     
#   }
#   
#   #查看p_val和OR值
#   IMR_p_val <- vector(mode = "numeric",length = 0)
#   IMR_odds_ratio <- vector(mode = "numeric",length = 0)
#   
#   for (i in 1:length(kkIMR)) {
#     IMR_p_val <- c(IMR_p_val,kkIMR[[i]]$p_val)
#     IMR_odds_ratio <- c(IMR_odds_ratio,kkIMR[[i]]$odds_ratio)
#   }
#   
#   length(IMR_p_val[which(IMR_p_val<0.01)])
#   IMR_odds_ratio <- IMR_odds_ratio[which(IMR_p_val<0.01)]
#   length(IMR_odds_ratio[which(IMR_odds_ratio<1)])
#   length(IMR_odds_ratio[which(IMR_odds_ratio>1)])
#   length(kkIMR)
#   #查看OR值大于10的gene
#   IMR_OR_greter_10 <- IMR_ID[which(IMR_odds_ratio>10)]
#   
# }
# 
# if(F){
#   ###GO_KEGG analysis
#   #
#   library(org.Hs.eg.db)
#   #gene <- DMR_OR_greter_10
#   gene <- IMR_OR_greter_10
#   allGene <- unique(human_region$id)
#   geneList <- bitr(allGene,fromType = "ENSEMBL",
#                    toType = "ENTREZID",
#                    OrgDb = org.Hs.eg.db)
#   gene_df <- bitr(gene, fromType = "ENSEMBL",
#                   toType = c("ENTREZID", "SYMBOL"),
#                   OrgDb = org.Hs.eg.db)
#   head(gene_df)
#   ego <- enrichGO(gene          = gene_df$ENTREZID,
#                   universe      = geneList$ENTREZID,
#                   OrgDb         = org.Hs.eg.db,
#                   ont           = "ALL",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.01,
#                   qvalueCutoff  = 0.05)
#   head(summary(ego))
#   dotplot(ego, showCategory=30)
# 
#   kk <- enrichKEGG(gene          = gene_df$ENTREZID,
#                    universe      = geneList$ENTREZID,
#                  organism = 'hsa',
#                  pvalueCutoff=0.05,
#                  pAdjustMethod="BH",
#                  qvalueCutoff=0.1)
#   head(kk)
#   dotplot(kk, showCategory=30)
# }
# 
# if(F){
#   #整体看DMR\IMR和CRE的富集分析
#   sum(width(GenomicRanges::intersect(h_dmr,hg38CRE,ignore.strand=T)))
#   sum(width(h_dmr))
#   sum(width(hg38CRE))
#   sum(width(hg38))
#   fisher.test(matrix(c(1287416,3683312,252960490,797919039),2,2))
#   
#   sum(width(GenomicRanges::intersect(h_imr,hg38CRE,ignore.strand=T)))
#   sum(width(h_imr))
#   sum(width(hg38CRE))
#   sum(width(hg38))
#   fisher.test(matrix(c(3493537,8694956,250754369,792907395),2,2))
#   
# }
# 
# if(F){
#   ##查看IMR所在基因与HouseKeeping基因的交集
#   library(readr)
#   housekeeoing <- read_delim("D:/R/R_project/housekeeoing.txt", 
#                              "\t", escape_double = FALSE, col_names = FALSE, 
#                              trim_ws = TRUE)
#   housekeeoing1 <- bitr(housekeeoing$X1,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db)
#   veen <- venn.diagram(list(IMR=unique(h_imr$id),HouseKeeping=housekeeoing1$ENSEMBL),
#                        filename = 'IMR_HouseKeeping.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue"),
#                        alpha = 0.5,
#                        label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        cat.dist = c(-0.05,-0.02)
#                        )
#   veen <- venn.diagram(list(IMR=unique(h_dmr$id),HouseKeeping=housekeeoing1$ENSEMBL),
#                        filename = 'DMR_HouseKeeping.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue"),
#                        alpha = 0.5,
#                        label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        cat.dist = c(-0.05,-0.02)
#   )
#   
# }
# 
# 
# #查看OR大于10的DMR长度分析
# DMR_OR_greter_10 <- which(DMR_odds_ratio>10)
# DMR1 <- GenomicRanges::GRanges()
# CRE1 <- GenomicRanges::GRanges()
# for (i in DMR_OR_greter_10) {
#   DMR1 <- c(DMR1,kkDMR[[i]]$DMR)
#   CRE1 <- c(CRE1,kkDMR[[i]]$CRE)
# }
# width(DMR1)
# 
# IMR_OR_greter_10 <- which(IMR_odds_ratio>10)
# IMR1 <- GenomicRanges::GRanges()
# CRE1 <- GenomicRanges::GRanges()
# for (i in IMR_OR_greter_10) {
#   IMR1 <- c(IMR1,kkIMR[[i]]$IMR)
#   CRE1 <- c(CRE1,kkIMR[[i]]$CRE)
# }
# width(IMR1)
# 
# 
# #查看OR值大于1和小于1的DMR的长度分布，有没有可能是长度过长导致的CRE富集
# OR_more1 <- EnsemblID_DMR[which(odds_ratio>1)]
# length(OR_more1)
# save(list = ls(),file = '../liver_ready.Rdata')
# rm(list=ls())
# #####################################################################################
# 
# 
# 

# 
# if(F){
# 
#   ###GO_KEGG analysis
#   #
#   # library(org.Hs.eg.db)
#   a = unique(human_region_brain[which(similarity_brain<0.2)]$id)
#   b = unique(human_region_brain[which(similarity_brain>0.9)]$id)
#   c <- intersect(a,b)
#   #gene <- setdiff(a,c)
#   gene <- setdiff(b,c)
#   #gene <- c
#   allGene <- unique(human_region_brain$id)
#   geneList <- bitr(allGene,fromType = "ENSEMBL",
#                    toType = "ENTREZID",
#                    OrgDb = org.Hs.eg.db)
#   gene_df <- bitr(gene, fromType = "ENSEMBL",
#                   toType = c("ENTREZID", "SYMBOL"),
#                   OrgDb = org.Hs.eg.db)
#   head(gene_df)
#   ego <- enrichGO(gene          = gene_df$ENTREZID,
#                   universe      = geneList$ENTREZID,
#                   OrgDb         = org.Hs.eg.db,
#                   ont           = "ALL",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.01,
#                   qvalueCutoff  = 0.05)
#   head(summary(ego))
#   dotplot(ego, showCategory=30)
# 
#   kk <- enrichKEGG(gene          = gene_df$ENTREZID,
#                    universe      = geneList$ENTREZID,
#                  organism = 'hsa',
#                  pvalueCutoff=0.05,
#                  pAdjustMethod="BH",
#                  qvalueCutoff=0.1)
#   head(kk)
#   dotplot(kk, showCategory=30)
# 
# }
# 
# if(F){
#   DMR_ID <- unique(h_dmr$id)
#   length(DMR_ID)
#   kkDMR <- list()
#   for (i in 1:length(DMR_ID)) {
#     #当前包含DMR基因
#     gene <- hg38[which(hg38$EnsemblID==DMR_ID[i])]
#     #当前基因的所有DMR
#     DMR <- GenomicRanges::intersect(h_dmr,gene,ignore.strand=TRUE)
#     #当前基因的所有CRE
#     CRE <- GenomicRanges::intersect(hg38CRE,gene,ignore.strand=TRUE)
#     a = ifelse(length(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))>0,sum(width(GenomicRanges::intersect(DMR,CRE,ignore.strand=TRUE))),0)
#     
#     b = sum(width(DMR)) - a 
#     
#     c = sum(width(CRE)) - a 
#     
#     d = width(gene) - (a+b+c)
#     
#     data <- data.frame(YesDMR=c(a,b),NoDMR=c(c,d))
#     rownames(data)=c("YesCRE","NoCRE")
#     fisher <- fisher.test(data)
#     kkDMR[[i]] <- list(data=data,p_val=fisher$p.value,odds_ratio=fisher$estimate,geneLength=gene,DMR=DMR,CRE=CRE)
#     
#   }
#   
#   #查看p_val和OR值
#   DMR_p_val <- vector(mode = "numeric",length = 0)
#   DMR_odds_ratio <- vector(mode = "numeric",length = 0)
#   
#   for (i in 1:length(kkDMR)) {
#     DMR_p_val <- c(DMR_p_val,kkDMR[[i]]$p_val)
#     DMR_odds_ratio <- c(DMR_odds_ratio,kkDMR[[i]]$odds_ratio)
#   }
#   
#   length(DMR_p_val[which(DMR_p_val<0.01)])
#   DMR_odds_ratio <- DMR_odds_ratio[which(DMR_p_val<0.01)]
#   length(DMR_odds_ratio[which(DMR_odds_ratio<1)])
#   length(DMR_odds_ratio[which(DMR_odds_ratio>10)])
#   length(kkDMR)
#   DMR_OR_greter_10 <- DMR_ID[which(DMR_odds_ratio>10)]
#   
# }
# if(F){
#   #查看liver的CRE与IMR、DMR的富集情况
#   library(readr)
#   ENCFF151YYY_ENCFF905FLR_7group <- read_delim("D:/R/R_project/ENCFF151YYY_ENCFF905FLR.7group.bed", 
#                                                "\t", escape_double = FALSE, col_names = FALSE, 
#                                                trim_ws = TRUE)
#   #去除Unclassified
#   ENCFF151YYY_ENCFF905FLR_7group <- ENCFF151YYY_ENCFF905FLR_7group[which(ENCFF151YYY_ENCFF905FLR_7group$X10!="Unclassified"),]
#   H3K4me3_H3K27ac <- GenomicRanges::GRanges(
#     seqnames = ENCFF151YYY_ENCFF905FLR_7group$X1,
#     ranges = IRanges::IRanges(start = ENCFF151YYY_ENCFF905FLR_7group$X2, end = ENCFF151YYY_ENCFF905FLR_7group$X3),
#     type = ENCFF151YYY_ENCFF905FLR_7group$X10
#   )
#   h_dmr = human_region[which(similarity<0.2)]
#   h_imr = human_region[which(similarity>0.9)]
#   #DMR与H3K4me3_H3K27ac富集
#   H3K4me3_H3K27ac <- GenomicRanges::intersect(H3K4me3_H3K27ac,hg38,ignore.strand=T)
#   h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3_H3K27ac,ignore.strand=T)))
#   b = sum(width(h_dmr)) - a
#   c = sum(width(H3K4me3_H3K27ac)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #IMR与H3K4me3_H3K27ac富集
#   h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3_H3K27ac,ignore.strand=T)))
#   b = sum(width(h_imr)) - a
#   c = sum(width(H3K4me3_H3K27ac)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
# 
#   
#   library(readr)
#   ENCFF641UEZ_7group <- read_delim("D:/R/R_project/ENCFF641UEZ.7group.bed", 
#                                    "\t", escape_double = FALSE, col_names = FALSE, 
#                                    trim_ws = TRUE)
#   
#   #去除Unclassified
#   ENCFF641UEZ_7group <- ENCFF641UEZ_7group[which(ENCFF641UEZ_7group$X10!="Unclassified"),]
#   
#   H3K4me3 <- GenomicRanges::GRanges(
#     seqnames = ENCFF641UEZ_7group$X1,
#     ranges = IRanges::IRanges(start = ENCFF641UEZ_7group$X2, end = ENCFF641UEZ_7group$X3),
#     type = ENCFF641UEZ_7group$X10
#   )
#   
#   #DMR与H3K4me3富集
#   H3K4me3 <- GenomicRanges::intersect(H3K4me3,hg38,ignore.strand=T)
#   h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3,ignore.strand=T)))
#   b = sum(width(h_dmr)) - a
#   c = sum(width(H3K4me3)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #IMR与H3K4me3富集
#   h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3,ignore.strand=T)))
#   b = sum(width(h_imr)) - a
#   c = sum(width(H3K4me3)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   library(readr)
#   ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1 <- read_delim("D:/R/R_project/ENCFF958VCA_ENCFF035NGT_ENCFF646FIY.7group (1).bed", 
#                                                               "\t", escape_double = FALSE, col_names = FALSE, 
#                                                               trim_ws = TRUE)
#   
#   
#   CTCF <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("CTCF",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#   H3K4me3 <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("H3K4me3",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#   H3K27ac <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("H3K27ac",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#   
#   H3K4me3 <- GenomicRanges::GRanges(
#     seqnames = H3K4me3$X1,
#     ranges = IRanges::IRanges(start = H3K4me3$X2, end = H3K4me3$X3),
#     type = H3K4me3$X10
#   )
#   
#   CTCF <- GenomicRanges::GRanges(
#     seqnames = CTCF$X1,
#     ranges = IRanges::IRanges(start = CTCF$X2, end = CTCF$X3),
#     type = CTCF$X10
#   )
#   H3K27ac <- GenomicRanges::GRanges(
#     seqnames = H3K27ac$X1,
#     ranges = IRanges::IRanges(start = H3K27ac$X2, end = H3K27ac$X3),
#     type = H3K27ac$X10
#   )
#   #DMR与H3K4me3富集
#   H3K4me3 <- GenomicRanges::intersect(H3K4me3,hg38,ignore.strand=T)
#   h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3,ignore.strand=T)))
#   b = sum(width(h_dmr)) - a
#   c = sum(width(H3K4me3)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #IMR与H3K4me3富集
#   h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3,ignore.strand=T)))
#   b = sum(width(h_imr)) - a
#   c = sum(width(H3K4me3)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #DMR与CTCF富集
#   CTCF <- GenomicRanges::intersect(CTCF,hg38,ignore.strand=T)
#   h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_dmr,CTCF,ignore.strand=T)))
#   b = sum(width(h_dmr)) - a
#   c = sum(width(CTCF)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #IMR与CTCF富集
#   h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_imr,CTCF,ignore.strand=T)))
#   b = sum(width(h_imr)) - a
#   c = sum(width(CTCF)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #DMR与H3K27ac富集
#   H3K27ac <- GenomicRanges::intersect(H3K27ac,hg38,ignore.strand=T)
#   h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_dmr,H3K27ac,ignore.strand=T)))
#   b = sum(width(h_dmr)) - a
#   c = sum(width(H3K27ac)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   #IMR与H3K27ac富集
#   h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#   a = sum(width(GenomicRanges::intersect(h_imr,H3K27ac,ignore.strand=T)))
#   b = sum(width(h_imr)) - a
#   c = sum(width(H3K27ac)) - a
#   d = sum(width(hg38)) -(a+b+c)
#   fisher.test(matrix(c(a,b,c,d),2,2))
#     
# }
# 
# if(F){
#   #查看HOX基因的甲基化模式
#   library(readr)
#   HOX <- read_delim("D:/R/R_project/HOX", "\t", 
#                     escape_double = FALSE, col_names = FALSE, 
#                     trim_ws = TRUE)
#   h_dmr = human_region[which(similarity<0.2)]
#   h_imr = human_region[which(similarity>0.9)]
#   
#   veen <- venn.diagram(list(IMR=unique(h_imr$id),HOX=HOX$X1),
#                        filename = 'DMR_HOX.tiff',
#                        height = 3000,
#                        width = 3000,
#                        resolution = 500,
#                        col = "transparent",
#                        fill = c("red","blue"),
#                        alpha = 0.5,
#                        label.col = c("darkred","white","darkblue"),
#                        cex = 2.5,
#                        fontfamily = "serif",
#                        fontface = "bold",
#                        cat.default.pos = "outer",#设置标签在圆外面
#                        cat.cex = 2,#外标签的字体大小
#                        cat.fontfamily = "serif",
#                        cat.dist = c(-0.05,-0.02)
#   )
#   
# }
# 
# 
# 
# ##查看liver的IMR所在的基因与brain的IMR所在的基因的交集
# save(h_imr_brain,h_dmr_brain,file = "../brain_IMR_DMR.Rdata")
# rm(list = ls())
# load('../brain_IMR_DMR.Rdata')
# load('../liver_ready.Rdata')
# 
# veen <- venn.diagram(list(liver=unique(h_imr$id),brain=unique(h_imr_brain$id)),
#                      filename = 'IMR_brain_liver.tiff',
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
#                      cat.dist = c(-0.05,-0.05)
# )
# 
# veen <- venn.diagram(list(liver=unique(h_dmr$id),brain=unique(h_dmr_brain$id)),
#                      filename = 'DMR_brain_liver.tiff',
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
#                      cat.dist = c(-0.05,-0.05)
# )
# 
# sum(width(GenomicRanges::intersect(h_imr,h_imr_brain,ignore.strand=T)))
# sum(width(h_imr))
# sum(width(h_imr_brain))
# 

