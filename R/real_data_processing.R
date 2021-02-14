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
# which(similarity>0.9)
# i = 6024
# plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#   obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
# plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#   obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))


#   load('../liver3k.Rdata')
#   ############长度统计
  a=width(human_region[which(similarity<0.2)])
  species = rep("Human DMR",NROW(a))
  a = data.frame(width=a,species=species)
  b=width(mouse_region[which(similarity<0.2)])
  species = rep("Mouse DMR",NROW(b))
  b = data.frame(width=b,species = species)
  c=width(human_region[which(similarity>0.9)])
  species = rep("Human CMR",NROW(c))
  c = data.frame(width=c,species=species)
  d=width(mouse_region[which(similarity>0.9)])
  species = rep("Mouse CMR",NROW(d))
  d = data.frame(width=d,species = species)
  e = rbind(a,b,c,d)

  ggplot(e, aes(x = factor(species), y = width, fill = factor(species))) +
    geom_boxplot(notch = TRUE,outlier.colour = NA) +
    scale_fill_brewer(palette = "Pastel2") +
    xlab(label = "Liver") +
    ylab(label = "Region Width") +
    labs(fill = "Categories") + 
    scale_y_continuous(limits = c(0,2000))


# if(F){
#   # #ChipSeeker Annotation
#   load(file = '../liver_30.Rdata')
#   h_dmr = human_region[which(similarity<0.2)]
#   m_dmr = mouse_region[which(similarity<0.2)]
#
#   h_imr = human_region[which(similarity>0.9)]
#   m_imr = mouse_region[which(similarity>0.9)]
#
#   covplot(h_dmr)
#   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#   h_dmr_anno = annotatePeak(h_dmr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
#   plotAnnoPie(h_dmr_anno)
#   vennpie(h_dmr_anno)
#   upsetplot(h_dmr_anno)
#
#   covplot(h_imr)
#   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#   h_imr_anno = annotatePeak(h_imr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
#   plotAnnoPie(h_imr_anno)
#   vennpie(h_imr_anno)
#   upsetplot(h_imr_anno)
# }
#
# if(F){
#   load('../liver_30_CHG.Rdata')
#   library(readr)
#   tissue_category_rna_liver_Tissue <- read_delim("D:/R/R_project/tissue_category_rna_liver_Tissue.tsv",
#                                                  "\t", escape_double = FALSE, trim_ws = TRUE)
#   View(tissue_category_rna_liver_Tissue)
#   liver_specific_gene <- tissue_category_rna_liver_Tissue$Ensembl
#   liver_specific_gene <- intersect(liver_specific_gene,unique(human_anno$id))
#   
#   #查看DMR或IMR所在的基因与肝脏组织特异性基因是否显著富集
#   liver_DMR_gene <- unique(human_region[which(similarity<0.2)]$id)
#   liver_IMR_gene <- unique(human_region[which(similarity>0.9)]$id)
#   liver_IMR_DMR_gene <- intersect(liver_IMR_gene,liver_DMR_gene)
#   liver_DMR_gene <- setdiff(liver_DMR_gene,liver_IMR_DMR_gene)
#   liver_IMR_gene <- setdiff(liver_IMR_gene,liver_IMR_DMR_gene)
#   save(liver_DMR_gene,liver_IMR_gene,file = "liver_DMR_IMR_gene.Rdata")
#   
#   veen <- venn.diagram(list(DMR_gene=liver_DMR_gene,specific_gene=liver_specific_gene),
#                        filename = 'liver_DMR_specific_gene.tiff',
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
#                        cat.dist = c(-0.04,-0.01)
#   )
#   veen <- venn.diagram(list(IMR=liver_IMR_gene,specific=liver_specific_gene),
#                        filename = 'liver_IMR_specific_gene.tiff',
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
#                        cat.dist = c(-0.04,-0.01)
#   )
#   a=length(intersect(liver_specific_gene,liver_IMR_gene))
#   b=length(liver_IMR_gene)-a
#   c=length(liver_specific_gene)-a
#   d=length(human_anno$id)-a-b-c
#   fisher.test(matrix(c(a,b,c,d),2,2),alternative = "less")
#   
#   a=length(intersect(liver_specific_gene,liver_DMR_gene))
#   b=length(liver_DMR_gene)-a
#   c=length(liver_specific_gene)-a
#   d=length(human_anno$id)-a-b-c
#   fisher.test(matrix(c(a,b,c,d),2,2),alternative = "greater")
#   
# }
#
# if(F){
#   #motif分析
#   load('../liver_3k.Rdata')
#   dmr <- human_region[which(similarity<0.2)]
#   id <- paste0(dmr$id,dmr$center)
#   chr <- as.vector(dmr@seqnames)
#   start <- start(dmr)
#   end <- end(dmr)
#   strand <- as.vector(strand(dmr))
#   DMR <- data.frame(id=id,chr=chr,start=start,end=end,strand=strand)
#   write.table(DMR,file = "../liver_DMR.txt",sep = "\t",row.names = F,col.names = F)
#   ##sed -i 's/"//g'
#   ##findMotifsGenome.pl liver_DMR.txt hg38 liver_DMR -size 200 -mask
#
#   imr <- human_region[which(similarity>0.9)]
#   id <- paste0(imr$id,imr$center)
#   chr <- as.vector(imr@seqnames)
#   start <- start(imr)
#   end <- end(imr)
#   strand <- as.vector(strand(imr))
#   IMR <- data.frame(id=id,chr=chr,start=start,end=end,strand=strand)
#   write.table(IMR,file = "../liver_IMR.txt",sep = "\t",row.names = F,col.names = F)
#   ##sed -i 's/"//g'
#   #findMotifsGenome.pl liver_IMR.txt hg38 liver_IMR -size 200 -mask
#
# }
#
# if(F){
#   #人类样本中不同组织的IMR和DMR
#   load('../human_liver_brain.Rdata')
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
#   )
#   veen <- venn.diagram(list(DMR=unique(h_dmr$id),HouseKeeping=housekeeoing1$ENSEMBL),
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
#   liver_DMR_gene <- unique(human_region[which(similarity<0.2)]$id)
#   liver_IMR_gene <- unique(human_region[which(similarity>0.9)]$id)
#   liver_IMR_DMR_gene <- intersect(liver_IMR_gene,liver_DMR_gene)
#   liver_DMR_gene <- setdiff(liver_DMR_gene,liver_IMR_DMR_gene)
#   liver_IMR_gene <- setdiff(liver_IMR_gene,liver_IMR_DMR_gene)
#   a=length(intersect(housekeeoing1$ENSEMBL,liver_IMR_gene))
#   b=length(liver_IMR_gene)-a
#   c=length(housekeeoing1$ENSEMBL)-a
#   d=length(human_anno$id)-a-b-c
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   a=length(intersect(housekeeoing1$ENSEMBL,liver_DMR_gene))
#   b=length(liver_DMR_gene)-a
#   c=length(housekeeoing1$ENSEMBL)-a
#   d=length(human_anno$id)-a-b-c
#   fisher.test(matrix(c(a,b,c,d),2,2))
# }
#
# #跨组织做富集分析，肝脏的IMR与脑的IMR做富集分析，肝脏的DMR与脑的DMR做富集分析
# if(F){
#   load('liver_DMR_IMR_gene.Rdata')
#   load('brain_DMR_IMR_gene.Rdata')
#   liver_brain_IMR_gene <- intersect(liver_IMR_gene,brain_IMR_gene)
#   houskeeping_gene <- housekeeoing1$ENSEMBL
#   a=length(intersect(houskeeping_gene,liver_brain_IMR_gene))
#   b=length(liver_brain_IMR_gene)-a
#   c=length(houskeeping_gene)-a
#   d=12837-a-b-c
#   fisher.test(matrix(c(a,b,c,d),2,2))
#   
#   
#   
# }
#

##查看DMR和IMR所在的基因交集情况
  # a = unique(human_region[which(similarity<0.2)]$id)
  # b = unique(human_region[which(similarity>0.9)]$id)
  # veen <- venn.diagram(list(DMR=a,CMR=b),
  #                      filename = 'liver_DMR_IMR_gene.tiff',
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
  #                      cat.dist = c(-0.03,-0.1)
  # )

# if(F){
#     ###GO_KEGG analysis
#     #
#     library(org.Hs.eg.db)
#     library(clusterProfiler)
#     a = unique(human_region[which(similarity<0.2)]$id)
#     b = unique(human_region[which(similarity>0.9)]$id)
#     c <- intersect(a,b)
#     gene <- setdiff(a,c)
#     #gene <- setdiff(b,c)
#     #gene <- c
#     allGene <- unique(human_region$id)
#     geneList <- bitr(allGene,fromType = "ENSEMBL",
#                      toType = "ENTREZID",
#                      OrgDb = org.Hs.eg.db)
#     gene_df <- bitr(gene, fromType = "ENSEMBL",
#                     toType = c("ENTREZID", "SYMBOL"),
#                     OrgDb = org.Hs.eg.db)
#     head(gene_df)
#     ego <- enrichGO(gene          = gene_df$ENTREZID,
#                     universe      = geneList$ENTREZID,
#                     OrgDb         = org.Hs.eg.db,
#                     # ont           = "ALL",
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05)
#     head(summary(ego))
#     dotplot(ego, showCategory=30)
#     
#     kk <- enrichKEGG(gene          = gene_df$ENTREZID,
#                      universe      = geneList$ENTREZID,
#                      organism = 'hsa',
#                      pvalueCutoff=0.05,
#                      pAdjustMethod="BH",
#                      qvalueCutoff=0.1)
#     head(kk)
#     dotplot(kk, showCategory=30)
# }

#if(F){
#     
#     #  #ChipSeeker Annotation
#     #   load(file = '../liver_3k.Rdata')
#     h_dmr = human_region[which(similarity<0.2)]
#     m_dmr = mouse_region[which(similarity<0.2)]
#
#     h_imr = human_region[which(similarity>0.9)]
#     m_imr = mouse_region[which(similarity>0.9)]
#
#     covplot(h_dmr)
#     txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#     h_dmr_anno = annotatePeak(h_dmr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
#     plotAnnoPie(h_dmr_anno)
#     vennpie(h_dmr_anno)
#     upsetplot(h_dmr_anno)
#
#     covplot(h_imr)
#     txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#     h_imr_anno = annotatePeak(h_imr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
#     plotAnnoPie(h_imr_anno)
#     vennpie(h_imr_anno)
#     upsetplot(h_imr_anno)
#     
#}

#       ############查看GTEx对应组织的top100表达基因
#       library(readr)
#       liverGTEx_Portal <- read_delim("D:/downloads/liverGTEx Portal.txt",
#                                      "\t", escape_double = FALSE, trim_ws = TRUE)
#       #substr(liverGTEx_Portal$`Gencode Id`,1,15)
#       veen <- venn.diagram(list(DMR=gene,TOP100=substr(liverGTEx_Portal$`Gencode Id`,1,15)),
#                            filename = 'liver_DMR_top100.tiff',
#                            height = 3000,
#                            width = 3000,
#                            resolution = 500,
#                            col = "transparent",
#                            fill = c("red","blue"),
#                            alpha = 0.5,
#                            label.col = c("darkred","white","darkblue"),
#                            cex = 2.5,
#                            fontfamily = "serif",
#                            fontface = "bold",
#                            cat.default.pos = "outer",#设置标签在圆外面
#                            cat.cex = 2,#外标签的字体大小
#                            cat.fontfamily = "serif",
#                            cat.dist = c(-0.05,-0.02)
#       )
#
#       if(F){
#         #查看liver的CRE与IMR、DMR的富集情况
#         library(readr)
#         HUMAN_3K <- read_delim("D:/R/R_project/HUMAN_3K.gtf",
#                                "\t", escape_double = FALSE, col_names = FALSE,
#                                trim_ws = TRUE)
#         hg38 <- GenomicRanges::GRanges(
#           seqnames = paste0("chr",HUMAN_3K$X1),
#           ranges = IRanges::IRanges(start = HUMAN_3K$X4, end = HUMAN_3K$X5),
#           EnsemblID = HUMAN_3K$X2
#         )
#
#         ENCFF151YYY_ENCFF905FLR_7group <- read_delim("D:/R/R_project/ENCFF151YYY_ENCFF905FLR.7group.bed",
#                                                      "\t", escape_double = FALSE, col_names = FALSE,
#                                                      trim_ws = TRUE)
#         #去除Unclassified
#         ENCFF151YYY_ENCFF905FLR_7group <- ENCFF151YYY_ENCFF905FLR_7group[which(ENCFF151YYY_ENCFF905FLR_7group$X10!="Unclassified"),]
#         H3K4me3_H3K27ac <- GenomicRanges::GRanges(
#           seqnames = ENCFF151YYY_ENCFF905FLR_7group$X1,
#           ranges = IRanges::IRanges(start = ENCFF151YYY_ENCFF905FLR_7group$X2, end = ENCFF151YYY_ENCFF905FLR_7group$X3),
#           type = ENCFF151YYY_ENCFF905FLR_7group$X10
#         )
#         h_dmr = human_region[which(similarity<0.2)]
#         h_imr = human_region[which(similarity>0.9)]
#         #DMR与H3K4me3_H3K27ac富集
#         H3K4me3_H3K27ac <- GenomicRanges::intersect(H3K4me3_H3K27ac,hg38,ignore.strand=T)
#         h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3_H3K27ac,ignore.strand=T)))
#         b = sum(width(h_dmr)) - a
#         c = sum(width(H3K4me3_H3K27ac)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #IMR与H3K4me3_H3K27ac富集
#         h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3_H3K27ac,ignore.strand=T)))
#         b = sum(width(h_imr)) - a
#         c = sum(width(H3K4me3_H3K27ac)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#
#         library(readr)
#         ENCFF641UEZ_7group <- read_delim("D:/R/R_project/ENCFF641UEZ.7group.bed",
#                                          "\t", escape_double = FALSE, col_names = FALSE,
#                                          trim_ws = TRUE)
#
#         #去除Unclassified
#         ENCFF641UEZ_7group <- ENCFF641UEZ_7group[which(ENCFF641UEZ_7group$X10!="Unclassified"),]
#
#         H3K4me3 <- GenomicRanges::GRanges(
#           seqnames = ENCFF641UEZ_7group$X1,
#           ranges = IRanges::IRanges(start = ENCFF641UEZ_7group$X2, end = ENCFF641UEZ_7group$X3),
#           type = ENCFF641UEZ_7group$X10
#         )
#
#         #DMR与H3K4me3富集
#         H3K4me3 <- GenomicRanges::intersect(H3K4me3,hg38,ignore.strand=T)
#         h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3,ignore.strand=T)))
#         b = sum(width(h_dmr)) - a
#         c = sum(width(H3K4me3)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #IMR与H3K4me3富集
#         h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3,ignore.strand=T)))
#         b = sum(width(h_imr)) - a
#         c = sum(width(H3K4me3)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         library(readr)
#         ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1 <- read_delim("D:/R/R_project/ENCFF958VCA_ENCFF035NGT_ENCFF646FIY.7group (1).bed",
#                                                                     "\t", escape_double = FALSE, col_names = FALSE,
#                                                                     trim_ws = TRUE)
#
#
#         CTCF <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("CTCF",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#         H3K4me3 <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("H3K4me3",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#         H3K27ac <- ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1[which(grepl("H3K27ac",ENCFF958VCA_ENCFF035NGT_ENCFF646FIY_7group_1$X10)),]
#
#         H3K4me3 <- GenomicRanges::GRanges(
#           seqnames = H3K4me3$X1,
#           ranges = IRanges::IRanges(start = H3K4me3$X2, end = H3K4me3$X3),
#           type = H3K4me3$X10
#         )
#
#         CTCF <- GenomicRanges::GRanges(
#           seqnames = CTCF$X1,
#           ranges = IRanges::IRanges(start = CTCF$X2, end = CTCF$X3),
#           type = CTCF$X10
#         )
#         H3K27ac <- GenomicRanges::GRanges(
#           seqnames = H3K27ac$X1,
#           ranges = IRanges::IRanges(start = H3K27ac$X2, end = H3K27ac$X3),
#           type = H3K27ac$X10
#         )
#         #DMR与H3K4me3富集
#         H3K4me3 <- GenomicRanges::intersect(H3K4me3,hg38,ignore.strand=T)
#         h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_dmr,H3K4me3,ignore.strand=T)))
#         b = sum(width(h_dmr)) - a
#         c = sum(width(H3K4me3)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #IMR与H3K4me3富集
#         h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_imr,H3K4me3,ignore.strand=T)))
#         b = sum(width(h_imr)) - a
#         c = sum(width(H3K4me3)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #DMR与CTCF富集
#         CTCF <- GenomicRanges::intersect(CTCF,hg38,ignore.strand=T)
#         h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_dmr,CTCF,ignore.strand=T)))
#         b = sum(width(h_dmr)) - a
#         c = sum(width(CTCF)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #IMR与CTCF富集
#         h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_imr,CTCF,ignore.strand=T)))
#         b = sum(width(h_imr)) - a
#         c = sum(width(CTCF)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #DMR与H3K27ac富集
#         H3K27ac <- GenomicRanges::intersect(H3K27ac,hg38,ignore.strand=T)
#         h_dmr <- GenomicRanges::intersect(h_dmr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_dmr,H3K27ac,ignore.strand=T)))
#         b = sum(width(h_dmr)) - a
#         c = sum(width(H3K27ac)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#         #IMR与H3K27ac富集
#         h_imr <- GenomicRanges::intersect(h_imr,hg38,ignore.strand=T)
#         a = sum(width(GenomicRanges::intersect(h_imr,H3K27ac,ignore.strand=T)))
#         b = sum(width(h_imr)) - a
#         c = sum(width(H3K27ac)) - a
#         d = sum(width(hg38)) -(a+b+c)
#         fisher.test(matrix(c(a,b,c,d),2,2))
#
#       }
#





