edgeR_data <- list("Genotype_short" = NA,
                   "Genotype_long" = NA,
                   "Fast_WT" = NA,
                   "Fast_KO" = NA,
                   "Interaction"= NA)

for (i in 1:5){
    edgeR_data[[i]]<- openxlsx::read.xlsx(here::here("data/edgeR.xlsx"),i) |> as.data.table()
}

edgeR_sig <- edgeR_data
for (i in 1:5){
    edgeR_sig[[i]]<-edgeR_sig[[i]] %>%
        dplyr::filter(FDR < 0.05)
    edgeR_sig[[i]]<-edgeR_sig[[i]]$SYMBOL
}
cpm_matrix <- openxlsx::read.xlsx(here("data/CPM_matrix.xlsx"))
metadata <- load_metadata("metadata.xlsx")
metadata <- metadata %>%
    dplyr::mutate(Group = paste(Genotype, Fast, sep = "_")) |>
    dplyr::arrange(Sample)

#####Upsetplot#####
order_upset <- c("Genotype_short", "Genotype_long", "Fast_WT","Fast_KO","Interaction")
upset_data <- data.frame("Gene"= edgeR_sig[[1]],
                         "Group"= names(edgeR_sig[1]))
for (i in 2:length(edgeR_sig)){
    upset_data <- dplyr::add_row(upset_data,"Gene"=edgeR_sig[[i]],
                         "Group"=names(edgeR_sig)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Gene))

rownames(upset_wide)<-upset_wide$Gene

upset_wide <- dplyr::select(upset_wide,-Gene)

upsetRNA <- ComplexUpset::upset(upset_wide, order_upset, name = "", sort_sets = F, themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=theme(text = element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA + ggplot2::ggtitle("Upsetplot")+ggplot2::theme(plot.title = element_text(size = 18,
                                                                                              hjust = 0.5,
                                                                                              vjust = 0.95))

upsetRNA
 # tiff("UpsetRNA.tif", unit = "cm", height = 25, width = 50, res = 600)
 # upsetRNA
 # grid::grid.text("Dif. expressed genes", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
 # dev.off()

#####code for GO analysis####

####GO on signficant genes####
go_sig_genes <- goAnalysis(edgeR_data, "MF", "")
#printGOterms(go_sig_genes)

go_sig_genes_cc <- goAnalysis(edgeR_data, "CC", "")
#printGOterms(go_sig_genes_cc)

for (i in 1:length(go_sig_genes)){
    go_sig_genes[[i]] <- clusterProfiler::setReadable(go_sig_genes[[i]],
                                                    OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
}

for (i in 1:length(go_sig_genes_cc)){
    go_sig_genes_cc[[i]] <- clusterProfiler::setReadable(go_sig_genes_cc[[i]],
                                                      OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
}

#extract data from the organelle inner membrane GO term for heatmap
go_sig_genes_cc$Genotype_long$Description[4]
organelle_inner_membrane <- go_sig_genes_cc$Genotype_long$geneID[4]
organelle_inner_membrane <- unlist(str_split(organelle_inner_membrane, "/"))

res_ox <- cpm_matrix %>%
    dplyr::filter(SYMBOL %in% organelle_inner_membrane) %>%
    dplyr::distinct(SYMBOL, .keep_all = T)
rownames(res_ox) <- res_ox$SYMBOL
res_ox <- res_ox %>%
    dplyr::select(-SYMBOL)

setup_ordered <- metadata
setup_ordered <- setup_ordered %>%
    dplyr::mutate(
        Group = dplyr::case_when(
            Group == "WT_18h" ~ "WT 18h",
            Group == "KO_18h" ~ "HNKO 18h",
            Group == "WT_5h" ~ "WT 5h",
            Group == "KO_5h" ~ "HNKO 5h"
        )
    )
order <- c("WT 5h", "HNKO 5h", "WT 18h", "HNKO 18h")

setup_ordered <- setup_ordered %>%
    dplyr::arrange(desc(Fast),desc(Genotype))

setup_ordered$Sample <- as.character(setup_ordered$Sample)
res_ox <- res_ox %>%
    dplyr::select(setup_ordered$Sample)

key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(Group)
rownames(key) <- setup_ordered$Sample
key$Group <- factor(key$Group, c("WT 5h", "HNKO 5h", "WT 18h", "HNKO 18h"))



Heatmap <- pheatmap::pheatmap(res_ox,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            scale = "row",
                            legend = T,
                            na_col = "white",
                            Colv = NA,
                            na.rm = T,
                            cluster_cols = F,
                            fontsize_row = 3,
                            fontsize_col = 8,
                            cellwidth = 7,
                            cellheight = 3,
                            annotation_col = key,
                            show_colnames = F,
                            show_rownames = T,
                            main = "Organelle Inner Membrane"
)

#####run GO analysis on data unique for various categories####
#preprocess and subset upset data
upset_data <- data.frame("Gene"= edgeR_sig[[1]],
                         "Group"= names(edgeR_sig[1]))
for (i in 2:length(edgeR_sig)){
    upset_data <- dplyr::add_row(upset_data,"Gene"=edgeR_sig[[i]],
                                 "Group"=names(edgeR_sig)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
long_genotype <- upset_data |>
    dplyr::filter(Group=="Genotype_long")
short_genotype <- upset_data |>
    dplyr::filter(Group=="Genotype_short")
fast_wt <- upset_data |>
    dplyr::filter(Group == "Fast_WT")
fast_ko <- upset_data |>
    dplyr::filter(Group=="Fast_KO")

#Create bg list for GO analysis
bg <- cpm_matrix$SYMBOL
bg_list <- clusterProfiler::bitr(
    bg,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

#analyze common genotype data

overlap <- long_genotype |>
    dplyr::filter(Gene %in% short_genotype$Gene)

genes <- unlist(str_split(overlap$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                       universe = bg_list$ENTREZID,
                                       OrgDb = org.Mm.eg.db,
                                       ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                         ont = "CC")
#####create heatmap on 100 random genes####
res_ox <- cpm_matrix %>%
    dplyr::filter(SYMBOL %in% genes) %>%
    dplyr::distinct(SYMBOL, .keep_all = T) |>
    dplyr::slice(1:100)
res_ox <- res_ox |>
    dplyr::filter(!is.na(SYMBOL))
rownames(res_ox) <- res_ox$SYMBOL
res_ox <- res_ox %>%
    dplyr::select(-SYMBOL)

setup_ordered <- metadata
setup_ordered <- setup_ordered %>%
    dplyr::mutate(
        Group = dplyr::case_when(
            Group == "WT_18h" ~ "WT 18h",
            Group == "KO_18h" ~ "HNKO 18h",
            Group == "WT_5h" ~ "WT 5h",
            Group == "KO_5h" ~ "HNKO 5h"
        )
    )
order <- c("WT 5h", "HNKO 5h", "WT 18h", "HNKO 18h")

setup_ordered <- setup_ordered %>%
    dplyr::arrange(desc(Fast),desc(Genotype))

setup_ordered$Sample <- as.character(setup_ordered$Sample)
res_ox <- res_ox %>%
    dplyr::select(setup_ordered$Sample)

key <- as.data.frame(setup_ordered)

key <- key %>%
    dplyr::select(Group)
rownames(key) <- setup_ordered$Sample
key$Group <- factor(key$Group, c("WT 5h", "HNKO 5h", "WT 18h", "HNKO 18h"))

Heatmap <- pheatmap::pheatmap(res_ox,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 3,
                              fontsize_col = 8,
                              cellwidth = 7,
                              cellheight = 3,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = "100 genes w. genotype regardless of fast"
)
######analysis of genes that's only different in long genotype####

unique_long <- long_genotype |>
    dplyr::filter(!Gene %in% short_genotype$Gene)

genes <- unlist(str_split(unique_long$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC")

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")



 # tiff(here::here("Data/figures/Heatmap_inner_membrane_long_fast.tif"), unit = "cm", height = 25, width = 50, res = 600)
 # Heatmap
 # dev.off()

#####Terms unique for long fast genotype####
unique_long <- long_genotype |>
    dplyr::filter(!Gene %in% short_genotype$Gene & !Gene %in% fast_wt$Gene & !Gene %in% fast_ko$Gene)
genes <- unlist(str_split(unique_long$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC")

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

organelle_inner_membrane <- heatmap_generator(goResults_CC$geneID[4], cpm_matrix, metadata, goResults_CC$Description[4])

  # tiff(here::here("Data/figures/Heatmap_inner_mito_membrane_long_fast_unique.tif"), unit = "cm", height = 25, width = 50, res = 600)
  # organelle_inner_membrane
  # dev.off()

#####Terms that change for long fast genotype and genotype effect with fast####
candidates <- long_genotype |>
    dplyr::filter(!Gene %in% short_genotype$Gene & !Gene %in% fast_wt$Gene & Gene %in% fast_ko$Gene)
genes <- unlist(str_split(candidates$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC")

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

genes_coerced <- paste(genes, collapse = "/")

long_genotype_fast_KO_overlap <- heatmap_generator(genes_coerced, cpm_matrix, metadata, "Candidates for dif. behaviour in HNKO")

  # tiff(here::here("Data/figures/Heatmap_long_genotype_fast_KO_overlap.tif"), unit = "cm", height = 25, width = 50, res = 600)
  # long_genotype_fast_KO_overlap
  # dev.off()

#####Genes that change in all conditions####
candidates <- long_genotype |>
    dplyr::filter(Gene %in% short_genotype$Gene & Gene %in% fast_wt$Gene & Gene %in% fast_ko$Gene)
genes <- unlist(str_split(candidates$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC")

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
goResults_MF <- clusterProfiler::setReadable(goResults_MF, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
genes_coerced <- paste(genes, collapse = "/")

all_changing_genes <- heatmap_generator_clustered(genes_coerced, cpm_matrix, metadata, "Genes behaving differently at all conditions")

  # tiff(here::here("Data/figures/changed_at_all_conditions.tif"), unit = "cm", height = 30, width = 50, res = 600)
  # all_changing_genes
  # dev.off()

#####Genes changed by fasting####
candidates <- fast_wt |>
    dplyr::filter(Gene %in% fast_ko$Gene)
genes <- unlist(str_split(candidates$Gene, "/"))

gene_list <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

goResults_MF <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_CC <- clusterProfiler::enrichGO(gene = gene_list$ENTREZID,
                                          universe = bg_list$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC")

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
goResults_MF <- clusterProfiler::setReadable(goResults_MF, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

 organelle_inner_membrane <- heatmap_generator(goResults_CC$geneID[1], cpm_matrix, metadata, goResults_CC$Description[1])
    # tiff(here::here("Data/figures/fasting_response.tif"), unit = "cm", height = 30, width = 50, res = 600)
    # organelle_inner_membrane
    # dev.off()

 #re-create MDS plots
 cpm_mds <- cpm_matrix |>
     dplyr::select(-ENSEMBL, -SYMBOL) |>
     dplyr::select(as.character(metadata$Sample))

 mdsData <- limma::plotMDS(as.matrix(cpm_mds), plot = F)
 mdsData <-
     mdsData$eigen.vectors %>% as.data.table() %>%
     dplyr::mutate(ID = metadata$Sample) %>%
     dplyr::mutate(Group = metadata$Group) %>%
     dplyr::select(ID, Group, V1, V2, V3)
 setnames(mdsData,
          c("V1", "V2", "V3", "ID", "Group"),
          c("dim1", "dim2", "dim3", "ID", "Group"))
 pBase <-
     ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
     geom_point(size = 5) +
     #geom_label(show.legend = FALSE, size = 5) +
     theme_bw()
pBase

#####Test reactomePA####

reactome_test <- ReactomePA::enrichPathway()

reactome_test <- reactomeAnalysis(edgeR_data, "")

edgeR_sig_entrez <- edgeR_sig
for (i in 1:length(edgeR_sig_entrez)){
    edgeR_sig_entrez[[i]] <- clusterProfiler::bitr(
        edgeR_sig_entrez[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[,2]

}

bg <- clusterProfiler::bitr(
    edgeR_data[[1]]$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

clusterTest <- clusterProfiler::compareCluster(geneClusters = edgeR_sig_entrez[-5],
                                fun = "enrichGO",
                                    universe = bg$ENTREZID,
                                OrgDb = org.Mm.eg.db,
                                ont = "MF")
clusterProfiler::dotplot(clusterTest, showCategory = 20)

#  tiff(here::here("Data/figures/clustercomparison.tif"), unit = "cm", height = 30, width = 50, res = 600)
# clusterTest
# dev.off()

