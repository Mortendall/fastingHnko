####data import and trimming ####


counts <- count_matrix_assembly("count_matrix.xlsx")

metadata <- load_metadata("metadata.xlsx")
metadata <- metadata %>% dplyr::mutate(Group = paste(Genotype, Fast, sep = "_"))
data_collection <- openxlsx::read.xlsx(here::here("data-raw/additionalData.xlsx"))
super_meta <- dplyr::left_join(metadata,data_collection, by = c("Sample"="ID")) |>
    dplyr::select(-"Group.y") |>
    dplyr::rename(Group = "Group.x") |>
    dplyr::mutate(Necrosis_grouped = dplyr::case_when(
        Necrosis >=1 & Genotype == "KO" ~ paste0("Necrosis",Group),
        Necrosis <1 & Genotype == "KO"~ paste0("Normal",Group),
        TRUE ~ Group
    ))
#openxlsx::write.xlsx(super_meta, here::here("data/super_meta.xlsx"))
design <- stats::model.matrix( ~0+Necrosis_grouped, super_meta)
colnames(design) <-
    stringr::str_remove_all(colnames(design), "\\(|\\)|Necrosis_grouped|:")



ctrsts <- limma::makeContrasts(
    Genotype_normal5h = NormalKO_5h - WT_5h,
    Genotype_necrosis5h = NecrosisKO_5h - WT_5h,
    NecrosisHNKO5h = NecrosisKO_5h - NormalKO_5h,
    FastingWT = WT_18h - WT_5h,
    FastingNormal = NormalKO_18h - NormalKO_5h,
    FastingNecrosis = NecrosisKO_18h - NecrosisKO_5h,
    Genotype_normal18h = NormalKO_18h - WT_18h,
    Genotype_necrosis18h = NecrosisKO_18h - WT_18h,
    NecrosisHNKO18h = NecrosisKO_18h - NormalKO_18h,
    levels = design)

super_meta$Sample <- as.character(super_meta$Sample)
counts <- counts %>% dplyr::select(super_meta$Sample)
all(metadata$Sample == colnames(counts))


group <- as.matrix(super_meta$Necrosis_grouped)
RNAseq <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(RNAseq, design = design, min.count = 20)
RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
RNAseq <- edgeR::calcNormFactors(RNAseq)
key <- clusterProfiler::bitr(rownames(RNAseq$counts), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
RNAseq <- edgeR::estimateDisp(RNAseq,design)
efit <- edgeR::glmQLFit(RNAseq, design)
dgeResults <- apply(ctrsts, 2, . %>%
                        edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                        edgeR::topTags(n = Inf, p.value = 1) %>%
                        magrittr::extract2("table") %>%
                        data.table::as.data.table(keep.rownames = TRUE))
dgeResults_annotated <- dgeResults
for (i in 1:length(dgeResults_annotated)){
    data.table::setnames(dgeResults_annotated[[i]],names(dgeResults_annotated[[i]])[1], "ENSEMBL")
    ens2symbol <-
        clusterProfiler::bitr(dgeResults_annotated[[i]]$ENSEMBL,
                              fromType = 'ENSEMBL',
                              toType = 'SYMBOL',
                              OrgDb = "org.Mm.eg.db")
    ens2symbol %<>% data.table %>% data.table::setkey(ENSEMBL)
    dgeResults_annotated[[i]] <- dplyr::full_join(dgeResults_annotated[[i]], ens2symbol)
}
#write.xlsx(dgeResults_annotated, file = here("data/edgeR_necrosisSubset.xlsx"), asTable = TRUE)
#####Make some additional analysis####
sheet_names <- openxlsx::getSheetNames(here::here("data/edgeR_necrosisSubset.xlsx"))
edgeR_data_sub <- vector(mode = "list", length = length(sheet_names))

for (i in 1:length(sheet_names)){
    edgeR_data_sub[[i]]<- openxlsx::read.xlsx(here::here("data/edgeR_necrosisSubset.xlsx"),i) |> as.data.table()
}
names(edgeR_data_sub)<-sheet_names
edgeR_sig_sub <- edgeR_data_sub
for (i in 1:length(sheet_names)){
    edgeR_sig_sub[[i]]<-edgeR_sig_sub[[i]] |>
        dplyr::filter(FDR < 0.05)
    edgeR_sig_sub[[i]]<-edgeR_sig_sub[[i]]$SYMBOL
}
cpm_matrix <- openxlsx::read.xlsx(here::here("data/CPM_matrix.xlsx"))
metadata <- openxlsx::read.xlsx(here::here("data/super_meta.xlsx"))

go_sig_genes_sub <- goAnalysis(edgeR_data_sub, "MF", "")
#printGOterms(go_sig_genes)

go_sig_genes_cc_sub <- goAnalysis(edgeR_data_sub, "CC", "")
go_sig_genes_bp_sub<-goAnalysis(edgeR_data_sub, "BP", "")
go_sig_genes_bp_sub_up<-goAnalysis(edgeR_data_sub, "BP", "Upregulated")
go_sig_genes_bp_sub_down<-goAnalysis(edgeR_data_sub, "BP", "Downregulated")

edgeR_sig_entrez <- edgeR_sig_sub
for (i in 1:length(edgeR_sig_entrez)){
    edgeR_sig_entrez[[i]] <- clusterProfiler::bitr(
        edgeR_sig_entrez[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[,2]

}

bg <-  clusterProfiler::bitr(
    edgeR_data_sub[[1]]$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)[,2]
go_compare_cc <- clusterProfiler::compareCluster(geneClusters = edgeR_sig_entrez,
                                fun = "enrichGO",
                                universe = bg,
                                OrgDb = org.Mm.eg.db,
                                ont = "CC",
                                readable = T)


 # tiff(filename = here::here("data/figures/go_cc_plot.tiff"), width = 20, height = 20, units = "cm", res = 200)
 # enrichplot::dotplot(go_compare_cc,
 #                     size= "count")+
 #     ggplot2::theme(axis.text.x = element_text(angle = 45,hjust = 1))+
 #    ggplot2::ggtitle("Cluster Comparison - Cellular Component")
 # dev.off()

go_compare_mf <- clusterProfiler::compareCluster(geneClusters = edgeR_sig_entrez,
                                                 fun = "enrichGO",
                                                 universe = bg,
                                                 OrgDb = org.Mm.eg.db,
                                                 ont = "MF",
                                                 readable = T)

 # tiff(filename = here::here("data/figures/go_mf_plot.tiff"), width = 20, height = 30, units = "cm", res = 200)
 # enrichplot::dotplot(go_compare_mf,
 #                     size= "count",
 #                     showCategory=4)+
 #     ggplot2::theme(axis.text.x = element_text(angle = 45,hjust = 1))+
 #     ggplot2::ggtitle("Cluster Comparison - Molecular Function (ShowCategory = 4)")
 # dev.off()

heatmap_test <- heatmap_generator_clustered_sub(go_sig_genes_bp_sub[["Genotype_necrosis5h"]]@result$geneID[11], cpm_matrix, metadata, go_sig_genes_bp_sub[["Genotype_necrosis5h"]]@result$Description[11])
tiff(filename = here::here("data/figures/heatmap_FAO_sub.tiff"), width = 20, height = 30, units = "cm", res = 200)
    heatmap_test
dev.off()



# printGOterms(go_sig_genes_bp_sub, "GOBP_all")
# printGOterms(go_sig_genes_bp_sub_up, "GOBP_up")
# printGOterms(go_sig_genes_bp_sub_down, "GOBP_down")
#printGOterms(go_sig_genes_cc_sub, "GOsub_all")
#printGOterms(go_sig_genes_sub, "GOsubMF_all")

#upsetplot
#####Upsetplot#####
order_upset <- names(edgeR_data_sub)
upset_data <- data.frame("Gene"= edgeR_sig_sub[[1]],
                         "Group"= names(edgeR_sig_sub[1]))
for (i in 2:length(edgeR_sig_sub)){
    upset_data <- dplyr::add_row(upset_data,"Gene"=edgeR_sig_sub[[i]],
                                 "Group"=names(edgeR_sig_sub)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Gene))

rowname_vector<-upset_wide$Gene

upset_wide <- dplyr::select(upset_wide,-Gene)
rownames(upset_wide)<-rowname_vector
order_upset <- order_upset[-9]
upsetRNA <- ComplexUpset::upset(upset_wide,
                                order_upset,
                                name = "",
                                sort_sets = "ascending",
                                min_size = 10,
                                themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=theme(text = element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA + ggplot2::ggtitle("")+ggplot2::theme(plot.title = element_text(size = 18,
                                                                                              hjust = 0.5,
                                                                                              vjust = 0.95))
#
# upsetRNA
#  tiff("UpsetRNA_subanalysis.tif", unit = "cm", height = 25, width = 50, res = 300)
#  upsetRNA
# grid::grid.text("Dif. expressed genes", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
#  dev.off()
