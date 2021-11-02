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
metadata <- metadata %>% dplyr::mutate(Group = paste(Genotype, Fast, sep = "_"))

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


