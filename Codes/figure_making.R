### figure making###
pacman::p_load("bayesbio")
library(bayesbio)
# BiocManager::install("pRolocdata")
library(tidyverse)
library(openxlsx)
library("pRolocdata")
data(hyperLOPITU2OS2018)
annot_hyperlopit <- hyperLOPITU2OS2018@featureData@data %>% rownames_to_column("Uniprot")
write_tsv(annot_hyperlopit, here::here("Datasets","Processed","annot_hyperlopit.tsv"))
read_tsv(here::here("Datasets","Processed","annot_hyperlopit.tsv"))
Jaccard_index_list<-function(list_of_vectors, max_jacc = 0.5,steps = 1){
    
    #removes top step most similar with other pathways and most similar until max_jacc is reached
    #samples as cols
    # list_of_vectors = terms_to_reduce
    row_max = 1
    # top_similar = 1
    # max_jacc = 0.3
    # steps = 1
    while(row_max>max_jacc){
        # list_of_vectors <- enrichment_df_BP %>% subset(`p.adjust`<0.05 ) %>% 
        #   dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
        #   pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
        index<-map(.x = list_of_vectors, ~.x %>% list(.) %>% rep(.,length(list_of_vectors)) %>% 
                       map2_dbl(.x = .,.y = list_of_vectors,~bayesbio::jaccardSets(.x,.y))) %>% 
            imap_dfr(.x = ., ~set_names(.x,names(list_of_vectors)) %>% enframe(name = "Pathway2",value = "JaccIndex") %>% mutate(Pathway1 = .y)) %>% 
            subset(JaccIndex != 1) %>% as.data.table() 
        setkey(index,Pathway1,Pathway2)
        index <- index[Pathway1>Pathway2]
        row_max <- index$JaccIndex %>% max()
        index <- index[JaccIndex == row_max][,top_similar:= fifelse(str_detect(Pathway1,"nucleo|eactive"),Pathway2,Pathway1)]
        
        top_similar <- index$top_similar
        # pivot_wider(names_from = "Pathway2", values_from = "JaccIndex") %>% 
        # column_to_rownames("Pathway1") %>% as.matrix()
        # diag(index) <- 0
        # row_max <- index %>% matrixStats::rowMaxs() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% .[1]
        #top_similar <- index %>% matrixStats::rowSums2() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% names() %>% .[1:steps]
        # top_similar <- c(top_similar)
        list_of_vectors[which(names(list_of_vectors)%in%top_similar)]<-NULL
        print(row_max)
    }
    list_of_vectors
}


txt_folder <-  here::here("Datasets","Raw",
                          "EC_118",
                          "txt")

columns_to_keep = c("Majority protein IDs",
                    "LFQ intensity","Potential","Reverse")

Samples_df <- data.frame(order = 1:12 %>% as.character() %>% paste0(" ",.),
                         raw = c("P11906","P11907","P11908","P11909","P11910","P11911","P11912","P11913","P11914","P11915","P11916","P11917"),
                         sample = c(paste0("N_UC_",1:3),
                                    paste0("N_C_",1:3),
                                    paste0("C_C_",1:3),
                                    paste0("C_UC_",1:3)))


Data_matrix <- fread(here::here(txt_folder,
                                "proteinGroups.txt")) %>%
    dplyr::select(matches(paste(columns_to_keep, collapse = "|"))) %>% 
    set_names(.,str_replace_all(colnames(.),Samples_df %>% pull(sample,order) %>% rev())) %>% 
    janitor::clean_names() %>% 
    subset(!((reverse == "+") | (potential_contaminant == "+"))) %>% 
    dplyr::select(-c(reverse,potential_contaminant)) %>% 
    mutate(majority_protein_i_ds  = str_remove_all(majority_protein_i_ds,";[:graph:]*")) %>% 
    column_to_rownames("majority_protein_i_ds") %>% 
    mutate(across(everything(),~na_if(.x,0))) %>% 
    set_names(.,str_remove_all(colnames(.),"lfq_intensity_")) %>% 
    glimpse()
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% 
    set_names(c("Uniprot","Type","ID"))

DEP_data <- DEP_DDA(Data_matrix,"Rito_Fabio_EC_118")
Data_matrices <- map(.x = 1:2,
                     ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118matrices.xlsx"), sheet = .x, rowNames = T)) %>% 
    set_names(getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118matrices.xlsx")))  

enrichment <- Data_matrices$Unimputted %>% dplyr::select(matches("uc")) %>% 
    rownames_to_column("Uniprot") %>% 
    pivot_longer(-Uniprot,names_to = "samples",values_to = "abundance") %>% 
    mutate(samples = str_remove_all(samples,"_.$")) %>% 
    group_by(Uniprot,samples) %>% 
    summarise(mean_abundance = mean(abundance, na.rm = T)) %>% 
    pivot_wider(names_from = "samples", values_from = mean_abundance) %>% na.omit() %>% 
    mutate(nuclear_enrichemnt = n_uc-c_uc) %>% 
    dplyr::select(Uniprot, nuclear_enrichemnt)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
install.packages('ggridges')
# clusterprofiler
enrichment_GSEA <- clusterProfiler::gseGO(geneList     = enrichment %>% pull(nuclear_enrichemnt, Uniprot) %>% sort(decreasing = T) ,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              keyType = "UNIPROT",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

    enrichplot::ridgeplot(enrichment_GSEA)+scale_fill_gradientn("p.adjust", colours = c( "#0073C2FF","#EFC000FF")) +ggtitle(glue::glue("n_uc_vs_c_uc CC-GSEA"))
    ggsave(here::here("Output","EC_118" , "n_uc_vs_c_uc CC-GSEA.pdf"), height = 20, width  = 15)

    enrichment %>% left_join(annot_hyperlopit %>% dplyr::select(Uniprot, `final.assignment`)) %>% na.omit() %>% 
        group_by(`final.assignment`) %>% 
        add_count() %>% subset(n>30) %>% subset(`final.assignment` != "unknown") %>%
        ggplot(aes(x = `final.assignment`, y = nuclear_enrichemnt))+
        geom_boxplot()+theme_bw()+
        ggtitle("Nuclear enrichment protocol enriches for nuclear proteins")
    ggsave(here::here("Output","EC_118" , "n_uc_vs_c_uc change.pdf"), height = 20, width  = 15)
    
Volcano_DFs <- map(.x = 1:2,
                   ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx"), sheet = .x)) |> 
    set_names(getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx")))
Volcano_DFs$n_c_vs_n_uc_diff |> 
    left_join(annot_hyperlopit) |> 
    ggplot(aes(x = `final.assignment`,y = log2_FC))+
    geom_boxplot()+geom_point(alpha = 0.3)+theme_bw()+
    ggtitle('n_c_vs_n_uc_diff')
ggsave(here::here("Output","n_c_vs_n_uc_diff_compartments.pdf"),width = 15, height =  12)
Volcano_DFs$c_c_vs_c_uc_diff |> 
    left_join(annot_hyperlopit) |> 
    ggplot(aes(x = `final.assignment`,y = log2_FC))+
    geom_boxplot()+geom_point(alpha = 0.3)+theme_bw()+
    ggtitle("c_c_vs_c_uc_diff")
ggsave(here::here("Output","c_c_vs_c_uc_diff_compartments.pdf"),width = 15, height =  12)

Significant_proteins <- c(Volcano_DFs$c_c_vs_c_uc_diff %>% 
                              subset(p.adj<0.05) %>% 
                              pull(Uniprot),
                          Volcano_DFs$n_c_vs_n_uc_diff %>% 
                              subset(p.adj<0.05) %>% 
                              pull(Uniprot) ) %>% unique()

ego2 <- gseGO(geneList     = Volcano_DFs$c_c_vs_c_uc_diff %>% pull(log2_FC, Uniprot) %>% sort(decreasing = T) ,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              keyType = "UNIPROT",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
if((ego2@result %>% nrow)>0){
    enrichplot::ridgeplot(ego2)+ ggtitle(glue::glue("c_c_vs_c_uc_diff ALL-GSEA"))
    ggsave(here::here("Output","EC_118" , glue::glue("c_c_vs_c_uc_diff","DDA_118", "ALL-GSEA.png")), height = 20, width  = 15)}


ego3 <- gseGO(geneList     = Volcano_DFs$n_c_vs_n_uc_diff %>% pull(log2_FC, Uniprot) %>% sort(decreasing = T) ,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              keyType = "UNIPROT",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
if((ego3@result %>% nrow)>0){
    enrichplot::ridgeplot(ego3, showCategory = 68)+
        ggtitle(glue::glue("n_c_vs_n_uc_diff"," ALL-GSEA"))
    ggsave(here::here("Output","EC_118", glue::glue("n_c_vs_n_uc_diff","DDA_118", "ALL-GSEA.png")), height = 20, width  = 15)}

terms_to_reduce <- ego3@result %>% subset(`p.adjust`<0.05 ) %>% 
    dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
    pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.3)

write.xlsx(ego3@result %>% subset(Description %in% names(reduced_terms)), here::here("Datasets","Processed","GSEA_ALL_reduced_terms.xlsx")) 
reduced_terms<- reduced_terms %>% purrr::imap_dfr(.x = .,~enframe(.x,name = "pathway", value = "Uniprot") %>% mutate(pathway = .y))
Volcano_DFs$n_c_vs_n_uc_diff %>% 
    left_join(reduced_terms %>% subset(str_detect(pathway,'ATP|cellular|matrix'))) %>% 
    mutate(significant = if_else(`p.adj`< 0.05 & abs(log2_FC)>1, T,F),
           Significant_alpha = if_else(significant ==T|!is.na(pathway),T,F)
           ) %>% 
    # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
    ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = pathway, alpha = Significant_alpha))+
    
    geom_point(data = . %>% subset(is.na(pathway)), size = 3.5)+
    geom_point(data = . %>% subset(!is.na(pathway)), size = 3.5)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(pathway)),size = 8)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(pathway)), size = 8)+
    
    
    scale_colour_manual(values = 
                            c("mitochondrial matrix" = "#6EA376",
                              "cellular amino acid metabolic process" = "#A95049" ,
                              # "bounding membrane of organelle" = "green",
                              "ATP metabolic process" ="#FFA100"),
    )+
    labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
    # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
    annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Confined","Unconfined")),size = 7, colour = "#8F9191")+
    lims(x = c(-3.19,3.19))+ theme_bw()+
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    #       text = element_text(size = 25),
    #       panel.background = element_blank()
    # )+
    ggtitle(glue::glue("Nuclear Changes Confinement"))
ggsave(here::here("Output","Annotated_GO_Proteomic_Volcano_nuclear.pdf"), height = 10, width = 10)
Volcano_DFs$c_c_vs_c_uc_diff %>% 
    left_join(reduced_terms %>% subset(str_detect(pathway,'ATP|cellular|matrix'))) %>% 
    mutate(significant = if_else(`p.adj`< 0.05 & abs(log2_FC)>1, T,F),
           Significant_alpha = if_else(significant ==T|!is.na(pathway),T,F)
    ) %>% 
    # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
    ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = pathway, alpha = Significant_alpha))+
    
    geom_point(data = . %>% subset(is.na(pathway)), size = 3.5)+
    geom_point(data = . %>% subset(!is.na(pathway)), size = 3.5)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(pathway)),size = 8)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(pathway)), size = 8)+
    
    
    scale_colour_manual(values = 
                            c("mitochondrial matrix" = "#6EA376",
                              "cellular amino acid metabolic process" = "#A95049" ,
                              # "bounding membrane of organelle" = "green",
                              "ATP metabolic process" ="#FFA100"),
    )+
    labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
    # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
    annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Confined","Unconfined")),size = 7, colour = "#8F9191")+
    lims(x = c(-3.19,3.19))+ theme_bw()+
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    #       text = element_text(size = 25),
    #       panel.background = element_blank()
    # )+
    ggtitle(glue::glue("Cytosol Changes Confinement"))
ggsave(here::here("Output","Annotated_GO_Proteomic_Volcano_cyto.pdf"), height = 10, width = 10)


# 
# ##### CORREP#####
# Correp_correlation <- correp(Data_matrices$Imputted,Significant_proteins)
# Correp_correlation$heatmap
# rowMeans <- Data_matrices$Imputted %>% rowMeans()
# 
# 
# 
# clusters_df <- Data_matrices$Imputted %>% 
#     mutate(across(everything(), ~.x-rowMeans)) %>% 
#     rownames_to_column("gene") %>% 
#     pivot_longer(-gene, names_to = "condition", values_to = "abundance") %>% 
#     inner_join(Correp_correlation$clusters %>% enframe("gene","cluster")) %>% 
#     mutate(cluster = as.character(cluster))
# 
# 
# 
# clusters_df_avg <- clusters_df %>% 
#     group_by(condition,cluster) %>% 
#     summarise(abundance = median(abundance))
# 
# level_order <- c("n_uc_1",  "n_uc_2", "n_uc_3",
#                  "n_c_1",  "n_c_2",  "n_c_3",
#                  "c_uc_1", "c_uc_2", "c_uc_3",
#                  "c_c_1",  "c_c_2",  "c_c_3")
# 
# ggplot(data=clusters_df,aes(x = factor(condition, level = level_order),  y= abundance, colour = cluster, group = gene))+
#     geom_line(alpha = 0.1)+
#     geom_point(alpha = 0.1)+
#     facet_wrap("cluster")+
#     geom_line(data=clusters_df_avg, aes(x=condition, y=abundance, group = cluster ), colour="grey10")+
#     ggtitle("Cluster Behaviour from CORREP") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     xlab("samples")
# ggsave(here::here("Output", "EC_118","Significant Protein Clusters reordered RG.png"))
# 
# 
# 
# list_of_diff_genes <- clusters_df %>% 
#     group_split(cluster) %>% 
#     set_names(.,map_chr(.x = .,~paste0("cluster_",.x %>% pull(cluster) %>% unique()))) %>% 
#     map(.x = .,~.x %>% pull(gene) %>% str_remove_all(";[:graph:]*$")%>% unique() )
# enrich_clusters <- function(gene_list){
#     enrichGO(gene         = gene_list,
#              OrgDb         = org.Hs.eg.db,
#              universe = rownames(Data_matrices$Imputted)%>% str_remove_all(";[:graph:]*$")%>% unique(),
#              keyType       = 'UNIPROT',
#              ont           = "BP",
#              pAdjustMethod = "BH",
#              pvalueCutoff  = 0.01,
#              qvalueCutoff  = 0.05) 
# }
# ck2 <- compareCluster(geneCluster = list_of_diff_genes, fun = enrich_clusters)
# dotplot(ck2)
# ck2@compareClusterResult %>% 
#     mutate(GeneRatio = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
#     group_by(Cluster) %>% 
#     arrange(p.adjust) %>% 
#     slice(1:20) %>% 
#     ggplot(aes(x = Cluster, y= Description, colour = p.adjust, size = GeneRatio))+
#     geom_point()+
#     ggtitle("CORREP Clusters Enrichment BP simplified")
# ggsave(here::here("Output", "EC_118","Significant Protein Clusters ORA.png"))
# ```



# RG modifications for plotting -------------------------------------------

### figure making###
pacman::p_load("bayesbio")
library(bayesbio)
# BiocManager::install("pRolocdata")
library(tidyverse)
library(openxlsx)
library("pRolocdata")
data(hyperLOPITU2OS2018)
annot_hyperlopit <- hyperLOPITU2OS2018@featureData@data %>% rownames_to_column("Uniprot")
Jaccard_index_list<-function(list_of_vectors, max_jacc = 0.5,steps = 1){
    
    #removes top step most similar with other pathways and most similar until max_jacc is reached
    #samples as cols
    # list_of_vectors = terms_to_reduce
    row_max = 1
    # top_similar = 1
    # max_jacc = 0.3
    # steps = 1
    while(row_max>max_jacc){
        # list_of_vectors <- enrichment_df_BP %>% subset(`p.adjust`<0.05 ) %>% 
        #   dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
        #   pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
        index<-map(.x = list_of_vectors, ~.x %>% list(.) %>% rep(.,length(list_of_vectors)) %>% 
                       map2_dbl(.x = .,.y = list_of_vectors,~bayesbio::jaccardSets(.x,.y))) %>% 
            imap_dfr(.x = ., ~set_names(.x,names(list_of_vectors)) %>% enframe(name = "Pathway2",value = "JaccIndex") %>% mutate(Pathway1 = .y)) %>% 
            subset(JaccIndex != 1) %>% as.data.table() 
        setkey(index,Pathway1,Pathway2)
        index <- index[Pathway1>Pathway2]
        row_max <- index$JaccIndex %>% max()
        index <- index[JaccIndex == row_max][,top_similar:= fifelse(str_detect(Pathway1,"nucleo|eactive"),Pathway2,Pathway1)]
        
        top_similar <- index$top_similar
        # pivot_wider(names_from = "Pathway2", values_from = "JaccIndex") %>% 
        # column_to_rownames("Pathway1") %>% as.matrix()
        # diag(index) <- 0
        # row_max <- index %>% matrixStats::rowMaxs() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% .[1]
        #top_similar <- index %>% matrixStats::rowSums2() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% names() %>% .[1:steps]
        # top_similar <- c(top_similar)
        list_of_vectors[which(names(list_of_vectors)%in%top_similar)]<-NULL
        print(row_max)
    }
    list_of_vectors
}


txt_folder <-  here::here("Datasets","Raw",
                          "EC_118",
                          "txt")

columns_to_keep = c("Majority protein IDs",
                    "LFQ intensity","Potential","Reverse")

Samples_df <- data.frame(order = 1:12 %>% as.character() %>% paste0(" ",.),
                         raw = c("P11906","P11907","P11908","P11909","P11910","P11911","P11912","P11913","P11914","P11915","P11916","P11917"),
                         sample = c(paste0("N_UC_",1:3),
                                    paste0("N_C_",1:3),
                                    paste0("C_C_",1:3),
                                    paste0("C_UC_",1:3)))


Data_matrix <- fread(here::here(txt_folder,
                                "proteinGroups.txt")) %>%
    dplyr::select(matches(paste(columns_to_keep, collapse = "|"))) %>% 
    set_names(.,str_replace_all(colnames(.),Samples_df %>% pull(sample,order) %>% rev())) %>% 
    janitor::clean_names() %>% 
    subset(!((reverse == "+") | (potential_contaminant == "+"))) %>% 
    dplyr::select(-c(reverse,potential_contaminant)) %>% 
    mutate(majority_protein_i_ds  = str_remove_all(majority_protein_i_ds,";[:graph:]*")) %>% 
    column_to_rownames("majority_protein_i_ds") %>% 
    mutate(across(everything(),~na_if(.x,0))) %>% 
    set_names(.,str_remove_all(colnames(.),"lfq_intensity_")) %>% 
    glimpse()
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% 
    set_names(c("Uniprot","Type","ID"))

DEP_data <- DEP_DDA(Data_matrix,"Rito_Fabio_EC_118")
Data_matrices <- map(.x = 1:2,
                     ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118matrices.xlsx"), sheet = .x, rowNames = T)) %>% 
    set_names(getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118matrices.xlsx")))  

enrichment <- Data_matrices$Unimputted %>% dplyr::select(matches("uc")) %>% 
    rownames_to_column("Uniprot") %>% 
    pivot_longer(-Uniprot,names_to = "samples",values_to = "abundance") %>% 
    mutate(samples = str_remove_all(samples,"_.$")) %>% 
    group_by(Uniprot,samples) %>% 
    summarise(mean_abundance = mean(abundance, na.rm = T)) %>% 
    pivot_wider(names_from = "samples", values_from = mean_abundance) %>% na.omit() %>% 
    mutate(nuclear_enrichemnt = n_uc-c_uc) %>% 
    dplyr::select(Uniprot, nuclear_enrichemnt)
enrichment_GSEA <- gseGO(geneList     = enrichment %>% pull(nuclear_enrichemnt, Uniprot) %>% sort(decreasing = T) ,
                         OrgDb        = org.Hs.eg.db,
                         ont          = "CC",
                         keyType = "UNIPROT",
                         #nPerm        = 1000,
                         minGSSize    = 50,
                         maxGSSize    = 500,
                         pvalueCutoff = 0.05,
                         verbose      = FALSE)

enrichplot::ridgeplot(enrichment_GSEA)+ggtitle(glue::glue("n_uc_vs_c_uc CC-GSEA"))
ggsave(here::here("Output","EC_118" , "n_uc_vs_c_uc CC-GSEA.pdf"), height = 20, width  = 15)

enrichment %>% left_join(annot_hyperlopit %>% dplyr::select(Uniprot, `final.assignment`)) %>% na.omit() %>% 
    group_by(`final.assignment`) %>% 
    add_count() %>% subset(n>30) %>% subset(`final.assignment` != "unknown") %>%
    ggplot(aes(x = `final.assignment`, y = nuclear_enrichemnt))+
    geom_boxplot()+theme_bw()+
    ggtitle("Nuclear enrichment protocol enriches for nuclear proteins")
ggsave(here::here("Output","EC_118" , "n_uc_vs_c_uc change.pdf"), height = 20, width  = 15)

Volcano_DFs <- map(.x = 1:2,
                   ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx"), sheet = .x)) |> 
    set_names(getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx")))
Volcano_DFs$n_c_vs_n_uc_diff |> 
    left_join(annot_hyperlopit) |> 
    ggplot(aes(x = `final.assignment`,y = log2_FC))+
    geom_boxplot()+geom_point(alpha = 0.3)+theme_bw()+
    ggtitle('n_c_vs_n_uc_diff')
ggsave(here::here("Output","n_c_vs_n_uc_diff_compartments.pdf"),width = 15, height =  12)
Volcano_DFs$c_c_vs_c_uc_diff |> 
    left_join(annot_hyperlopit) |> 
    ggplot(aes(x = `final.assignment`,y = log2_FC))+
    geom_boxplot()+geom_point(alpha = 0.3)+theme_bw()+
    ggtitle("c_c_vs_c_uc_diff")
ggsave(here::here("Output","c_c_vs_c_uc_diff_compartments.pdf"),width = 15, height =  12)

Significant_proteins <- c(Volcano_DFs$c_c_vs_c_uc_diff %>% 
                              subset(p.adj<0.05) %>% 
                              pull(Uniprot),
                          Volcano_DFs$n_c_vs_n_uc_diff %>% 
                              subset(p.adj<0.05) %>% 
                              pull(Uniprot) ) %>% unique()

ego2 <- gseGO(geneList     = Volcano_DFs$c_c_vs_c_uc_diff %>% pull(log2_FC, Uniprot) %>% sort(decreasing = T) ,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              keyType = "UNIPROT",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
if((ego2@result %>% nrow)>0){
    enrichplot::ridgeplot(ego2)+ ggtitle(glue::glue("c_c_vs_c_uc_diff ALL-GSEA"))
    ggsave(here::here("Output","EC_118" , glue::glue("c_c_vs_c_uc_diff","DDA_118", "ALL-GSEA.png")), height = 20, width  = 15)}


ego3 <- gseGO(geneList     = Volcano_DFs$n_c_vs_n_uc_diff %>% pull(log2_FC, Uniprot) %>% sort(decreasing = T) ,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              keyType = "UNIPROT",
              #nPerm        = 1000,
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
if((ego3@result %>% nrow)>0){
    enrichplot::ridgeplot(ego3, showCategory = 68)+
        ggtitle(glue::glue("n_c_vs_n_uc_diff"," ALL-GSEA"))
    ggsave(here::here("Output","EC_118", glue::glue("n_c_vs_n_uc_diff","DDA_118", "ALL-GSEA.png")), height = 20, width  = 15)}

terms_to_reduce <- ego3@result %>% subset(`p.adjust`<0.05 ) %>% 
    dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
    pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
reduced_terms <- Jaccard_index_list(terms_to_reduce,max_jacc = 0.3)

write.xlsx(ego3@result %>% subset(Description %in% names(reduced_terms)), here::here("Datasets","Processed","GSEA_ALL_reduced_terms.xlsx")) 
reduced_terms<- reduced_terms %>% purrr::imap_dfr(.x = .,~enframe(.x,name = "pathway", value = "Uniprot") %>% mutate(pathway = .y))
Volcano_DFs$n_c_vs_n_uc_diff %>% 
    left_join(reduced_terms %>% subset(str_detect(pathway,'ATP|cellular|matrix'))) %>% 
    mutate(significant = if_else(`p.adj`< 0.05 & abs(log2_FC)>1, T,F),
           Significant_alpha = if_else(significant ==T|!is.na(pathway),T,F)
    ) %>% 
    # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
    ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = pathway, alpha = Significant_alpha))+
    
    geom_point(data = . %>% subset(is.na(pathway)), size = 3.5)+
    geom_point(data = . %>% subset(!is.na(pathway)), size = 3.5)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(pathway)),size = 8)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(pathway)), size = 8)+
    
    
    scale_colour_manual(values = 
                            c("mitochondrial matrix" = "#6EA376",
                              "cellular amino acid metabolic process" = "#A95049" ,
                              # "bounding membrane of organelle" = "green",
                              "ATP metabolic process" ="#FFA100"),
    )+
    labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
    # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
    annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Confined","Unconfined")),size = 7, colour = "#8F9191")+
    lims(x = c(-3.19,3.19))+ theme_bw()+
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    #       text = element_text(size = 25),
    #       panel.background = element_blank()
    # )+
    ggtitle(glue::glue("Nuclear Changes Confinement"))
ggsave(here::here("Output","Annotated_GO_Proteomic_Volcano_nuclear.pdf"), height = 10, width = 10)
Volcano_DFs$c_c_vs_c_uc_diff %>% 
    left_join(reduced_terms %>% subset(str_detect(pathway,'ATP|cellular|matrix'))) %>% 
    mutate(significant = if_else(`p.adj`< 0.05 & abs(log2_FC)>1, T,F),
           Significant_alpha = if_else(significant ==T|!is.na(pathway),T,F)
    ) %>% 
    # dplyr::mutate(Metabolic_library = dplyr::if_else(Metabolic_library == "Non-Metabolic", "Non-Metabolic","Metabolic")) %>% 
    ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = pathway, alpha = Significant_alpha))+
    
    geom_point(data = . %>% subset(is.na(pathway)), size = 3.5)+
    geom_point(data = . %>% subset(!is.na(pathway)), size = 3.5)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  is.na(pathway)),size = 8)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T&  !is.na(pathway)), size = 8)+
    
    
    scale_colour_manual(values = 
                            c("mitochondrial matrix" = "#6EA376",
                              "cellular amino acid metabolic process" = "#A95049" ,
                              # "bounding membrane of organelle" = "green",
                              "ATP metabolic process" ="#FFA100"),
    )+
    labs(x = "Fold Change on Chromatin", y = "-log10(pvalue)")+
    # annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
    annotate("text", x = c(-1.5,1.5), y=0, label = rev(c("Confined","Unconfined")),size = 7, colour = "#8F9191")+
    lims(x = c(-3.19,3.19))+ theme_bw()+
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    #       text = element_text(size = 25),
    #       panel.background = element_blank()
    # )+
    ggtitle(glue::glue("Cytosol Changes Confinement"))
ggsave(here::here("Output","Annotated_GO_Proteomic_Volcano_cyto.pdf"), height = 10, width = 10)


# 
# ##### CORREP#####
# Correp_correlation <- correp(Data_matrices$Imputted,Significant_proteins)
# Correp_correlation$heatmap
# rowMeans <- Data_matrices$Imputted %>% rowMeans()
# 
# 
# 
# clusters_df <- Data_matrices$Imputted %>% 
#     mutate(across(everything(), ~.x-rowMeans)) %>% 
#     rownames_to_column("gene") %>% 
#     pivot_longer(-gene, names_to = "condition", values_to = "abundance") %>% 
#     inner_join(Correp_correlation$clusters %>% enframe("gene","cluster")) %>% 
#     mutate(cluster = as.character(cluster))
# 
# 
# 
# clusters_df_avg <- clusters_df %>% 
#     group_by(condition,cluster) %>% 
#     summarise(abundance = median(abundance))
# 
# level_order <- c("n_uc_1",  "n_uc_2", "n_uc_3",
#                  "n_c_1",  "n_c_2",  "n_c_3",
#                  "c_uc_1", "c_uc_2", "c_uc_3",
#                  "c_c_1",  "c_c_2",  "c_c_3")
# 
# ggplot(data=clusters_df,aes(x = factor(condition, level = level_order),  y= abundance, colour = cluster, group = gene))+
#     geom_line(alpha = 0.1)+
#     geom_point(alpha = 0.1)+
#     facet_wrap("cluster")+
#     geom_line(data=clusters_df_avg, aes(x=condition, y=abundance, group = cluster ), colour="grey10")+
#     ggtitle("Cluster Behaviour from CORREP") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     xlab("samples")
# ggsave(here::here("Output", "EC_118","Significant Protein Clusters reordered RG.png"))
# 
# 
# 
# list_of_diff_genes <- clusters_df %>% 
#     group_split(cluster) %>% 
#     set_names(.,map_chr(.x = .,~paste0("cluster_",.x %>% pull(cluster) %>% unique()))) %>% 
#     map(.x = .,~.x %>% pull(gene) %>% str_remove_all(";[:graph:]*$")%>% unique() )
# enrich_clusters <- function(gene_list){
#     enrichGO(gene         = gene_list,
#              OrgDb         = org.Hs.eg.db,
#              universe = rownames(Data_matrices$Imputted)%>% str_remove_all(";[:graph:]*$")%>% unique(),
#              keyType       = 'UNIPROT',
#              ont           = "BP",
#              pAdjustMethod = "BH",
#              pvalueCutoff  = 0.01,
#              qvalueCutoff  = 0.05) 
# }
# ck2 <- compareCluster(geneCluster = list_of_diff_genes, fun = enrich_clusters)
# dotplot(ck2)
# ck2@compareClusterResult %>% 
#     mutate(GeneRatio = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
#     group_by(Cluster) %>% 
#     arrange(p.adjust) %>% 
#     slice(1:20) %>% 
#     ggplot(aes(x = Cluster, y= Description, colour = p.adjust, size = GeneRatio))+
#     geom_point()+
#     ggtitle("CORREP Clusters Enrichment BP simplified")
# ggsave(here::here("Output", "EC_118","Significant Protein Clusters ORA.png"))
# ```

