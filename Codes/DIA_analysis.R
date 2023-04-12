#### Etop DIA ####
pacman::p_load( here, tidyverse, DEP,pheatmap,sva,imp4p,proDA)
               org.Hs.eg.db,clusterProfiler,ggridges,SubCellBarCode,eulerr,scales,data.table,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)
#Change the inference to Protein level
#Double run through NN
#### Loading####
source(here::here("Codes","functions.R"))
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
Human_hsa <- HUMAN_9606 %>% 
    subset(Type == "KEGG") %>% 
    mutate(ID = str_remove_all(ID,"hsa:")) %>% 
    dplyr::select(-Type) %>% 
    subset(!duplicated(Uniprot))

output_folder = here::here("Output","EC_118")
contaminant_uniprots <- seqinr::read.fasta(here::here(
    "Datasets","Raw","contaminants.fasta"),seqtype ="AA") %>% names() 
#####
list_files_to_analyse <- list(Rito_DIA = "Rito_trial_EC_118_report.tsv")

Samples_df <- data.frame(#order = 1:12 %>% as.character() %>% paste0(" ",.),
                         run = c("P11906","P11907","P11908","P11909","P11910","P11911","P11912","P11913","P11914","P11915","P11916","P11917"),
                         Condition = c(paste0("NUC_",1:3),
                                    paste0("NC_",1:3),
                                    paste0("CC_",1:3),
                                    paste0("CUC_",1:3))) 

Methods_DIA <- map(.x = list_files_to_analyse,
                   ~Load_DIA_NN_Data(.x,Samples_df))
Methods_DIA_DEP <- map2(.x = Methods_DIA ,
                        .y = list_files_to_analyse,
                        ~DEP_DIA(.x,.y))

#### DIA DDA correlation
Volcano_Df_DDA <-  map(.x = 1:2,
                       ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx"), sheet = .x, rowNames = F)) %>% 
    set_names(openxlsx::getSheetNames(here::here("Datasets","Processed", "Rito_Fabio_EC_118volcano_DFs.xlsx"))) 
Volcano_Df_DIA <-  map(.x = 1:2,
                       ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Rito_trial_DIA_EC_118_report.tsvvolcano_DFs.xlsx"), sheet = .x, rowNames = F)) %>% 
    set_names(openxlsx::getSheetNames(here::here("Datasets","Processed", "Rito_trial_DIA_EC_118_report.tsvvolcano_DFs.xlsx"))) 
Cyto <- list(DDA = Volcano_Df_DDA$c_c_vs_c_uc_diff ,
             DIA = Volcano_Df_DIA$cc_vs_cuc_diff) %>% 
    imap(.x = .,~.x %>%  subset(!duplicated(Single_Uniprot)) %>% 
             dplyr::select(Single_Uniprot,log2_FC,significant,Imputted_comparison) %>% remove_rownames %>% 
             column_to_rownames("Single_Uniprot") %>% set_names(paste0(colnames(.),.y)) %>%
                 rownames_to_column("Uniprot")) %>% purrr::reduce(inner_join, by = "Uniprot")
Cyto %>% 
    mutate(Significant = case_when(
        
        significantDDA == T & significantDIA == T ~ "Both",
        significantDDA == T | significantDIA == T ~ "One",
        TRUE ~ "None"
    )) %>%  
    mutate(Imputed = case_when(
        
        Imputted_comparisonDDA == T | Imputted_comparisonDIA == T ~ "Imputted",
        TRUE ~ "Non-Imputted"
    )) %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                  subset(!duplicated(Uniprot)) ) %>% 
    ggplot(aes(x = log2_FCDDA,y = log2_FCDIA, colour = Significant, label = ID, alpha= Imputed))+
    geom_point()+
    ggrepel::geom_label_repel(data = . %>% subset(Significant == "Both"))+
    scale_alpha_manual(values = c("Imputted" = 0.5,"Non-Imputted" = 1))
    xlim(-9,9)+ylim(-9,9)
ggsave(here::here("Output","EC_118","Cyto_DIA_DDA_agreement.png"))

Nucl <- list(DDA = Volcano_Df_DDA$n_c_vs_n_uc_diff ,
             DIA = Volcano_Df_DIA$nc_vs_nuc_diff) %>% 
    imap(.x = .,~.x %>%  subset(!duplicated(Single_Uniprot)) %>% 
             dplyr::select(Single_Uniprot,log2_FC,significant,Imputted_comparison) %>% remove_rownames %>% 
             column_to_rownames("Single_Uniprot") %>% set_names(paste0(colnames(.),.y)) %>%
             rownames_to_column("Uniprot")) %>% reduce(inner_join, by = "Uniprot")
Nucl %>% 
    mutate(Significant = case_when(
        
        significantDDA == T & significantDIA == T ~ "Both",
        significantDDA == T | significantDIA == T ~ "One",
        TRUE ~ "None"
    )) %>%  
    mutate(Imputed = case_when(
        
        Imputted_comparisonDDA == T | Imputted_comparisonDIA == T ~ "Imputted",
        TRUE ~ "Non-Imputted"
    )) %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                  subset(!duplicated(Uniprot)) ) %>% 
    ggplot(aes(x = log2_FCDDA,y = log2_FCDIA, colour = Significant, label = ID, alpha= Imputed))+
    geom_point()+
    ggrepel::geom_label_repel(data = . %>% subset(Significant == "Both"))+
    scale_alpha_manual(values = c("Imputted" = 0.5,"Non-Imputted" = 1))

ggsave(here::here("Output","EC_118","Nucl_DIA_DDA_agreement.png"))

