library(tidyverse)
library(openxlsx)
library(data.table)
library(clusterProfiler)
library("pRolocdata")
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(simplifyEnrichment)



# Nuclear enrichment of organelles plot -------------------------------------------
data(hyperLOPITU2OS2018)
annot_hyperlopit <- hyperLOPITU2OS2018@featureData@data %>% rownames_to_column("Uniprot")

Volcano_DFs <- map(.x = 1:2,
                   ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx"), sheet = .x)) |> 
    set_names(getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx")))

order=Volcano_DFs$n_c_vs_n_uc_diff %>%  
    left_join(annot_hyperlopit) %>%
    dplyr::filter(is.na(final.assignment)!=TRUE,
                  final.assignment!="unknown") %>%
    group_by(final.assignment) %>% 
    summarise(m.log2fc = median(log2_FC)) %>% 
    arrange(-m.log2fc) %>% 
    ungroup()

#Creating enrichment of organelles, showing Mitochondria as most upregulated
organelle_plot_df=Volcano_DFs$n_c_vs_n_uc_diff %>%  
    left_join(annot_hyperlopit) %>%
    dplyr::filter(is.na(final.assignment)!=TRUE,
                  final.assignment!="unknown") %>%
    mutate(across(final.assignment, factor, levels=rev(order$final.assignment)),
           type="Enriched (log2 FC > 2)",
           type=case_when(final.assignment!='MITOCHONDRION'~"Normal",
                          TRUE~as.character(type)))  %>% 
    group_by(final.assignment) %>% 
    mutate(med = median(log2_FC)) %>% 
    ungroup()

enrichment_organelle = ggplot(organelle_plot_df, aes(y = final.assignment, x = log2_FC))+
    geom_boxplot(aes(fill=med), outlier.shape=NA)+
    geom_point(shape=19,size=2,fill='black', alpha=0.3)+
    scale_fill_gradient(low="#dddddd", high="#0073C2FF")+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        # legend.position = "none",
        # legend.position = c(0.75,0.1),
        axis.title.y = element_blank())+
    xlim(-2, 5)

pdf("C:/Users/rghose/OneDrive - CRG - Centre de Regulacio Genomica/Rito_Fabio/Img_Analysis_Data/Organelle_enrichment.pdf", 
    width = 5, height = 7.5)
enrichment_organelle
dev.off()

# Circular barchart for process enrichment --------------------------------
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
    plot_ego3_df = ego3@result %>% dplyr::select(ONTOLOGY, Description, NES, p.adjust) %>%
        arrange(p.adjust)}

#Separating the results into GO --> MF, CC, BP and calculating semantic matrix based on terms
ego3@result %>% head(75) %>%  group_by(ONTOLOGY) %>% tally()

ego3_bp = ego3@result %>% head(75) %>% dplyr::filter(ONTOLOGY=='BP') #49terms
ego3_bp_semanticMatrix = GO_similarity(ego3_bp$ID)
Heatmap(ego3_bp_semanticMatrix)

ego3_mf = ego3@result %>% head(75) %>% dplyr::filter(ONTOLOGY=='MF') #6terms 
ego3_mf_semanticMatrix = GO_similarity(ego3_mf$ID)
Heatmap(ego3_mf_semanticMatrix)

ego3_cc = ego3@result %>% head(75) %>% dplyr::filter(ONTOLOGY=='CC') #22terms :: Shows main Mitochondrial terms clustered
ego3_cc_semanticMatrix = GO_similarity(ego3_cc$ID)
Heatmap(ego3_cc_semanticMatrix)

#Using heatmaps to check the semantic clustering. 

ht=draw(Heatmap(ego3_cc_semanticMatrix, row_km=4))
rht=row_order(ht)

ego3_bp_toADD = ego3_bp %>% mutate(toADD = case_when(grepl("electron",Description)==TRUE ~ "ADD",
                                                     grepl("mito",Description)==TRUE ~ "ADD",
                                                     TRUE~as.character(Description))) %>% 
    dplyr::filter(toADD == "ADD") %>% 
    dplyr::select(Description, NES, p.adjust) %>% 
    mutate(cluster=4)

#Using the row order from Heatmap to establish clusters.
ego3_cc_cluster = rbind(ego3_cc %>% 
                            dplyr::select(Description, NES, p.adjust) %>%
                            dplyr::filter(row_number() %in% rht[[1]]) %>% 
                            mutate(cluster=1),
                        ego3_cc %>% 
                            dplyr::select(Description, NES, p.adjust) %>%
                            dplyr::filter(row_number() %in% rht[[2]]) %>% 
                            mutate(cluster=2),
                        ego3_cc %>% 
                            dplyr::select(Description, NES, p.adjust) %>%
                            dplyr::filter(row_number() %in% rht[[3]]) %>% 
                            mutate(cluster=3),
                        ego3_cc %>% 
                            dplyr::select(Description, NES, p.adjust) %>%
                            dplyr::filter(row_number() %in% rht[[4]]) %>% 
                            mutate(cluster=4)
) %>% 
    rbind(data.frame(Description=c(" "),
                     NES=c(0),
                     p.adjust=c(0),
                     cluster=c(5))) %>% 
    dplyr::filter(Description!='mitochondrion') %>% 
    arrange(cluster) %>% 
    mutate(Description=factor(Description, levels=Description))

#For manually ordering, save and upload data.
write.csv(ego3_cc_cluster,
          "C:/Users/rghose/OneDrive - CRG - Centre de Regulacio Genomica/Rito_Fabio/Datasets/Processed/circ_bar_df.csv", 
          quote=F)

plot_df<-read.csv("C:/Users/rghose/OneDrive - CRG - Centre de Regulacio Genomica/Rito_Fabio/Datasets/Processed/circ_bar_df.csv") %>% 
    mutate(Description=gsub(" ", "\n", Description),
           Description=factor(Description, levels=Description),
           )
# x=reorder(str_wrap(ego3_cc_cluster$Description, 5), ego3_cc_cluster$NES)
circ_gg<-ggplot(plot_df %>% arrange(cluster), 
       aes(x=Description, y=NES, fill=p.adjust))+
    geom_rect(aes(xmin="mitochondrial\nprotein-containing\ncomplex",
                                  xmax="mitochondrial\nmembrane",
                                  ymin=0, ymax=max(NES)+0.1),
              fill="#ADD8E6")+
    # geom_hline(aes(yintercept = y), data.frame(y = c(0.5,1,1.5)), color = "grey") +
    geom_col(color="black", width=0.8) +
    geom_segment(aes(y = 0, xend = Description, yend = max(NES)),linetype = "dashed", color = "black") +   
    coord_polar()+
    theme_bw()+
    scale_y_continuous(limits = c(-0.5, 2.4),
                       # expand = c(0, 0),
                       breaks = c(0, 0.5, 1, 1.5, 2)) +
    scale_fill_gradientn("p.adjust", colours = c( "#0073C2FF","#EFC000FF")) +
    theme_minimal()+
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "gray12", size = 10),
        legend.position = "bottom")
    

pdf("C:/Users/rghose/OneDrive - CRG - Centre de Regulacio Genomica/Rito_Fabio/Img_Analysis_Data/Circ_GO_CC.pdf", 
    width = 7.5, height = 7.5)
circ_gg
dev.off()





# Volcano plot showing mitochondria and ATP enrichment --------------------
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
                            c("mitochondrial matrix" = "#0073C2FF",
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
