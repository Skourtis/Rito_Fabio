Load_DIA_NN_Data <- function(report_DIA_tsv_file, Samples_df){
     Samples_df = Samples_df
   report_DIA_tsv_file = list_files_to_analyse$Rito_DIA
    read_tsv(here::here("Datasets","Raw",report_DIA_tsv_file)) %>% 
        as.data.frame() %>% 
        subset(Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01 ) %>% janitor::clean_names() %>% 
        dplyr::select(matches("pg_max_lfq|run|protein_group|^genes$")) %>% 
        distinct() %>% 
        # mutate(old_pg =protein_group) %>%  
        subset(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))) %>% 
        remove_rownames() %>%
        mutate(run = run %>% str_match("-([:graph:]{6})-.$") %>% .[,2] %>% tolower()) %>%  
        left_join(Samples_df %>% 
                      mutate(across(everything(),tolower))) %>% 
        pivot_wider(-run, names_from = Condition ,values_from = pg_max_lfq) %>% 
        column_to_rownames("protein_group")%>% 
        dplyr::select(any_of(tolower(Samples_df$Condition))) %>% 
        as.matrix()}


DEP_DIA <- function(input_matrix,dataset_name){
    #this required a not normalised and not logged matrix with column names of condition in _rep format
    dataset_name = list_files_to_analyse$Rito_DIA
     input_matrix <- Methods_DIA$Rito_DIA
    
    # as.data.frame() 
    # mutate(across(everything(),~.x^2)) %>% 
    input_matrix <- input_matrix %>% as.data.frame()
    set.seed(2023)
    experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:graph:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                      subset(Type == "Gene_Name") %>% 
                      dplyr::select(-Type) %>% 
                      subset(!duplicated(Uniprot)), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
    if((experimental_design_DIA %>% group_by(condition) %>% 
        dplyr::summarise(Total_Rep = length(unique(replicate))) %>% 
        pull(Total_Rep) %>% min())>3){
        n_accepted_NA_per_condition  = 3
        ComBAT = T
    }else{
        n_accepted_NA_per_condition  = 1
        ComBAT = F
    }
    
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_overlap ",dataset_name,".pdf")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".pdf")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    png(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
    plot_missval(data_filt)
    dev.off()
    data_norm <- normalize_vsn(data_filt)
    data_norm@assays@data@listData[[1]] <-  DEqMS::equalMedianNormalization(data_norm@assays@data@listData[[1]])
    # data_norm@assays@data@listData[[1]] <-  proDA::median_normalization(data_norm@assays@data@listData[[1]])
    
    # data_norm@assays@data@listData[[1]] <- input_matrix
    
    DEP::meanSdPlot(data_norm)
    ggsave(here::here(output_folder,glue::glue("normalize_vsn ",dataset_name,".pdf")))
    
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    
    pca_res <- prcomp(data_norm@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=32) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="top") +
        ggtitle(dataset_name)+ 
        theme(plot.title = element_text(size = 20))
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.pdf")))
    
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
        data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)
        # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
        # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition))
        if(ComBAT == T){
            batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
            
            multiple_imputation <-  cbind(impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 2],#methodMNAR = "impute.pa",
                                                    conditions = experimental_design_DIA$condition[batch == 2] %>% as.factor(),
                                                    repbio = experimental_design_DIA$replicate[batch == 2] %>% as.factor()),
                                          impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 1],#methodMNAR = "impute.pa",
                                                    conditions = experimental_design_DIA$condition[batch == 1] %>% as.factor(),
                                                    repbio = experimental_design_DIA$replicate[batch == 1] %>% as.factor())) 
            
        }else{multiple_imputation <-  impute.mi(tab = data_norm@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                                conditions = experimental_design_DIA$condition %>% as.factor(),
                                                repbio = experimental_design_DIA$replicate %>% as.factor())
        
        
        }
        
        rownames(multiple_imputation) <- rownames(data_norm@assays@data@listData[[1]])
        colnames(multiple_imputation) <- colnames(data_norm@assays@data@listData[[1]])
        data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
        data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm@assays@data@listData[[1]]),]
    }else{
        data_imp <- data_norm
    }
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".pdf")))
    
    if(ComBAT == T){
        edata = data_imp@assays@data@listData[[1]]
        data_norm@assays@data@listData[[1]] <- data_norm@assays@data@listData[[1]] %>%
            as.data.frame() %>% 
            dplyr::select(!matches("dmso_(4|5)")) %>% 
            as.matrix()
        
        
        
        # parametric adjustment
        combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
        
        
        input_matrix_batch <-2^combat_edata1 %>% as.data.frame() %>% 
            dplyr::select(!matches("dmso_(4|5)"))
        experimental_design_DIA <-  data.frame(
            label = colnames(input_matrix_batch),
            condition =  str_remove_all(colnames(input_matrix_batch),"_[:graph:]*$"),
            replicate = str_remove_all(colnames(input_matrix_batch),"^[:graph:]*_") %>% as.numeric()
        )
        data_unique_Etop <- input_matrix_batch %>% rownames_to_column("name") %>% 
            left_join(HUMAN_9606 %>% 
                          subset(Type == "Gene_Name") %>% 
                          dplyr::select(-Type) %>% 
                          subset(!duplicated(Uniprot)), 
                      by = c("name" = "Uniprot")) %>%
            subset(!duplicated(name))
        
        Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix_batch))# get LFQ column numbers
        data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
        data_filt <- filter_missval(data_se, thr = 1)
        data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
        data_norm_batch <- normalize_vsn(data_filt)
        plot_normalization(data_se,data_norm_batch)
        data_norm_batch@assays@data@listData[[1]][is.na(data_norm@assays@data@listData[[1]])] <- NA
        
        
        if(data_norm_batch@assays@data@listData[[1]] %>% is.na() %>% any()){
            data_imp <- DEP::impute(data_norm_batch, fun = "MinProb", q = 0.01)
            # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
            # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition)
            batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
            
            multiple_imputation <-  impute.mi(tab = data_norm_batch@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                              conditions = experimental_design_DIA$condition %>% as.factor(),
                                              repbio = experimental_design_DIA$replicate %>% as.factor())
            
            
            rownames(multiple_imputation) <- rownames(data_norm_batch@assays@data@listData[[1]])
            colnames(multiple_imputation) <- colnames(data_norm_batch@assays@data@listData[[1]])
            data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
            data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm_batch@assays@data@listData[[1]]),]
        }else{
            data_imp <- data_norm_batch
        }
        pca_res <- prcomp(data_imp@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE) 
        #plot_missval(data_filt)
        # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
        
        var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
        pca_res$x %>% 
            as.data.frame %>%
            rownames_to_column("Sample") %>% 
            mutate(Condition = str_remove(Sample,"_.$")) %>% 
            ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
            ggrepel::geom_label_repel()+
            theme_bw(base_size=32) + 
            labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                 y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
            theme(legend.position="top") +
            ggtitle(glue::glue("Batch_Corrected",dataset_name))+ 
            theme(plot.title = element_text(size = 20))
        ggsave(here::here(output_folder,glue::glue(dataset_name,"batch PCA.pdf")))
        
        
        
    }
    data_diff_all_contrasts <- DEP::test_diff(data_imp,  type = "manual",
                                              test = c("cc_vs_cuc", "nc_vs_nuc"))
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0)
    
    #plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
    #ggsave(here::here(output_folder,glue::glue("Sample_correlation ",dataset_name,".pdf")))
    #plot_volcano(dep, contrast = "T0_vs_T24", label_size = 2, add_names = TRUE)
    # plot_heatmap(dep, type = "centered", kmeans = TRUE, 
    #              k = 6, col_limit = 4, show_row_names = FALSE,
    #              indicate = c("condition", "replicate"))
    # plot_single(dep, proteins = c("P05386","Q9H9B4"))
    
  
   
    Significant_protein_list <- data.frame(ProteinGroup = dep@elementMetadata$name,
                                           CYTO = dep@elementMetadata$cc_vs_cuc_significant,
                                           NUCLE = dep@elementMetadata$nc_vs_nuc_significant) %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot)))
    
    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        #mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        #left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #left_join(Interesting_proteins) %>% 
        #subset(!is.na(Behaviour)) %>% 
        #pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
        #group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        #mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        #dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>% dplyr::select(-Significant) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #   mutate(Condition = str_remove_all(Condition,"_.")) %>% 
        # ungroup() %>% 
        #   group_by(ProteinGroup,Condition) %>% 
        # dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = T)) %>% 
        pivot_wider(names_from = "Condition",values_from = "Abundance") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
        ungroup %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        subset(!is.na(ID)) %>% 
        column_to_rownames("ID") %>%
        dplyr::select(where(is.numeric)) %>%
        mutate(Rowmean = rowMeans(.),
               across(where(is.numeric),~.x- Rowmean)) %>% 
        dplyr::select(-Rowmean) %>% 
        as.matrix() 
    paletteLength <- 20
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins, na.rm = T)/paletteLength, max(Significant_proteins, na.rm = T), length.out=floor(paletteLength/2)))
    
    png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
    pheatmap::pheatmap(Significant_proteins[,experimental_design_DIA$label %>% sort()],cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", cluster_rows = T,
                       # scale = "row",
                       main = glue::glue(dataset_name, " \nSignificant Proteins- imputted"),color=myColor, breaks=myBreaks)
    dev.off()
    # if(names(dev.cur()) != "null device"){dev.off()}
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
    Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
         # i = (dep@elementMetadata %>% names() %>% str_subset("diff"))[1]
        contrast <- str_remove_all(i,"_diff")
        print(contrast)
        conditions <- contrast %>% str_match("([:graph:]*)_vs_([:graph:]*)$") %>% .[2:3]
        non_missing_in_all_comparison <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames(),
                                           data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames()) %>% unique()
        Imputted <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames(),
                      data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames()) %>% unique()
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Uniprot = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
            subset(Uniprot %in% non_missing_in_all_comparison) %>% 
            
            mutate(Imputted_comparison = Uniprot %in% Imputted,
                   Single_Uniprot = Uniprot %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>%
            left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type), by  = c("Single_Uniprot" = "Uniprot")) 
        volcano_df %>% ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, alpha = significant))+
            geom_point()+
            geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
            ggrepel::geom_label_repel(data = . %>% subset(significant == T& Imputted_comparison == F))+
            annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
            
            ggtitle(glue::glue("Diff Present on Chromatin", contrast),
                    subtitle = dataset_name)
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_significant",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
        
        Comparisons_list[[i]] <- volcano_df
        

        ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T) %>%
                        set_names(.,str_remove_all(names(.),";[:graph:]*$")),
                      OrgDb        = org.Hs.eg.db,
                      ont          = "MF",
                      keyType = "UNIPROT",
                      #nPerm        = 1000,
                      minGSSize    = 50,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
        if((ego3@result %>% nrow)>0){
            enrichplot::ridgeplot(ego3 %>% simplify(), showCategory = 68)+
                ggtitle(glue::glue(dataset_name," MF-GSEA",i))
            ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," MF-GSEA.pdf")), height = 20, width  = 15)}
        ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T)%>%
                        set_names(.,str_remove_all(names(.),";[:graph:]*$")),
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      keyType = "UNIPROT",
                      #nPerm        = 1000,
                      minGSSize    = 50,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
        if((ego3@result %>% nrow)>0){

            enrichplot::ridgeplot(ego3 %>% simplify(), showCategory = 68)+
                ggtitle(glue::glue(dataset_name,"BP-GSEA",i))
            ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," BP-GSEA.pdf")), height = 20, width  = 15)}

        gene_list <- volcano_df %>%
          dplyr::select(Single_Uniprot,log2_FC) %>%
          left_join(Human_hsa, by = c("Single_Uniprot" ="Uniprot")) %>%
          na.omit()%>%arrange(-log2_FC) %>%   pull(log2_FC,ID)

        kk2 <- gseKEGG(geneList     = gene_list,
                       organism     = 'hsa',
                       minGSSize    = 10,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
        if((kk2@result %>% nrow)>0){

          enrichplot::ridgeplot(kk2, showCategory = 100)+
            ggtitle(glue::glue(dataset_name," KEGG ",i))
          ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," KEGG.pdf")), height = 20, width  = 15)}


        mkk2 <- gseMKEGG(geneList = gene_list,
                         organism = 'hsa',
                         minGSSize = 10,
                         pvalueCutoff = 0.05)
        if((mkk2@result %>% nrow)>0){

          enrichplot::ridgeplot(mkk2, showCategory = 68)+
            ggtitle(glue::glue(dataset_name," MKEGG ",i))
          ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," MKEGG.pdf")), height = 20, width  = 15)}

        
    }
    
     data_matrices <- list(Imputted = data_imp@assays@data@listData[[1]],
                          Unimputted = data_norm@assays@data@listData[[1]], 
                          DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff")))
    write.xlsx(data_matrices[1:2], here::here("Datasets","Processed",glue::glue(dataset_name, "matrices.xlsx")), overwrite = T,rowNames = T)
    write.xlsx(data_matrices[[3]], here::here("Datasets","Processed",glue::glue(dataset_name, "volcano_DFs.xlsx")), overwrite = T)
    
}
Quality_Control <- function(x){
    ### takes as input the txt output folder of MaxQuant 
    ### in the form of the here::here function
    ### and produces a QC report
    pacman::p_temp("PTXQC")
    require(methods)
    r = createReport(x)
}
DEP_DDA <- function(input_matrix,dataset_name){
    #this required a not normalised and not logged matrix with column names of condition in _rep format
    
    
    
     dataset_name = "Rito_Fabio_EC_118"
     # 
     input_matrix <- Data_matrix
    
    # as.data.frame() 
    # mutate(across(everything(),~.x^2)) %>% 
    input_matrix <- input_matrix %>% as.data.frame()
    set.seed(2023)
    experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:digit:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                      subset(Type == "Gene_Name") %>% 
                      dplyr::select(-Type) %>% 
                      subset(!duplicated(Uniprot)), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
    if((experimental_design_DIA %>% group_by(condition) %>% 
        dplyr::summarise(Total_Rep = length(unique(replicate))) %>% 
        pull(Total_Rep) %>% min())>3){
        n_accepted_NA_per_condition  = 3
        ComBAT = T
    }else{
        n_accepted_NA_per_condition  = 1
        ComBAT = F
    }
    
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_overlap ",dataset_name,".pdf")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".pdf")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    pdf(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".pdf")), width = 10, height = 10) 
    plot_missval(data_filt)
    dev.off()
    data_norm <- normalize_vsn(data_filt)
    data_norm@assays@data@listData[[1]] <-  DEqMS::equalMedianNormalization(data_norm@assays@data@listData[[1]])
        
        # proDA::median_normalization(data_norm@assays@data@listData[[1]])
    
    # data_norm@assays@data@listData[[1]] <- input_matrix
    
    # DEP::meanSdPlot(data_norm)
    # ggsave(here::here(output_folder,glue::glue("normalize_vsn ",dataset_name,".pdf")))
    data_norm_ = data_norm
    data_norm_@assays@data@listData[[1]] = data_norm@assays@data@listData[[1]]+mean(
        data_se@assays@data@listData[[1]],na.rm = T)
    plot_normalization(data_se, data_norm_)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    
    pca_res <- prcomp(data_norm@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=32) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="top") +
        ggtitle(dataset_name)+ 
        theme(plot.title = element_text(size = 20))
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.pdf")))
    
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
        data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)
        # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
        # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition))
        if(ComBAT == T){
            batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
            
            multiple_imputation <-  cbind(impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 2],#methodMNAR = "impute.pa",
                                                    conditions = experimental_design_DIA$condition[batch == 2] %>% as.factor(),
                                                    repbio = experimental_design_DIA$replicate[batch == 2] %>% as.factor()),
                                          impute.mi(tab = data_norm@assays@data@listData[[1]][,batch == 1],#methodMNAR = "impute.pa",
                                                    conditions = experimental_design_DIA$condition[batch == 1] %>% as.factor(),
                                                    repbio = experimental_design_DIA$replicate[batch == 1] %>% as.factor())) 
            
        }else{
            multiple_imputation <-  impute.mi(tab = data_norm@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                                conditions = experimental_design_DIA$condition %>% as.factor(),
                                                repbio = experimental_design_DIA$replicate %>% as.factor())
        
        
        }
        
        rownames(multiple_imputation) <- rownames(data_norm@assays@data@listData[[1]])
        colnames(multiple_imputation) <- colnames(data_norm@assays@data@listData[[1]])
        data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
        data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm@assays@data@listData[[1]]),]
    }else{
        data_imp <- data_norm
    }
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".pdf")))
    
    if(ComBAT == T){
        edata = data_imp@assays@data@listData[[1]]
        data_norm@assays@data@listData[[1]] <- data_norm@assays@data@listData[[1]] %>%
            as.data.frame() %>% 
            dplyr::select(!matches("dmso_(4|5)")) %>% 
            as.matrix()
        
        
        
        # parametric adjustment
        combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
        
        
        input_matrix_batch <-10^combat_edata1 %>% as.data.frame() %>% 
            dplyr::select(!matches("dmso_(4|5)"))
        experimental_design_DIA <-  data.frame(
            label = colnames(input_matrix_batch),
            condition =  str_remove_all(colnames(input_matrix_batch),"_[:graph:]*$"),
            replicate = str_remove_all(colnames(input_matrix_batch),"^[:graph:]*_") %>% as.numeric()
        )
        data_unique_Etop <- input_matrix_batch %>% rownames_to_column("name") %>% 
            left_join(HUMAN_9606 %>% 
                          subset(Type == "Gene_Name") %>% 
                          dplyr::select(-Type) %>% 
                          subset(!duplicated(Uniprot)), 
                      by = c("name" = "Uniprot")) %>%
            subset(!duplicated(name))
        
        Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix_batch))# get LFQ column numbers
        data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
        data_filt <- filter_missval(data_se, thr = 1)
        data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
        data_norm_batch <- normalize_vsn(data_filt)
        plot_normalization(data_se,data_norm_batch)
        data_norm_batch@assays@data@listData[[1]][is.na(data_norm@assays@data@listData[[1]])] <- NA
        
        
        if(data_norm_batch@assays@data@listData[[1]] %>% is.na() %>% any()){
            data_imp <- DEP::impute(data_norm_batch, fun = "MinProb", q = 0.01)
            # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>n_accepted_NA_per_condition)
            # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<=(n_accepted_NA_per_condition)
            batch = if_else(str_match(colnames(input_matrix),"_(.)$")[,2] %>% as.numeric() <4,2,1)
            
            multiple_imputation <-  impute.mi(tab = data_norm_batch@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                              conditions = experimental_design_DIA$condition %>% as.factor(),
                                              repbio = experimental_design_DIA$replicate %>% as.factor())
            
            
            rownames(multiple_imputation) <- rownames(data_norm_batch@assays@data@listData[[1]])
            colnames(multiple_imputation) <- colnames(data_norm_batch@assays@data@listData[[1]])
            data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
            data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm_batch@assays@data@listData[[1]]),]
        }else{
            data_imp <- data_norm_batch
        }
        pca_res <- prcomp(data_imp@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE) 
        #plot_missval(data_filt)
        # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
        
        var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
        
        pca_res$x %>% 
            as.data.frame %>%
            rownames_to_column("Sample") %>% 
            mutate(Condition = str_remove(Sample,"_.$")) %>% 
            ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
            ggrepel::geom_label_repel()+
            labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                 y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
            theme(legend.position="top") +
            ggtitle(glue::glue("Batch_corrected ",dataset_name))+ 
            theme(plot.title = element_text(size = 20))
        ggsave(here::here(output_folder,glue::glue(dataset_name,"batch PCA.pdf")))
        
        
        
    }
    data_diff_all_contrasts <- DEP::test_diff(data_imp,  type = "manual",
                                              test = c("c_c_vs_c_uc", "n_c_vs_n_uc"))
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 1)
    
    #plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
    #ggsave(here::here(output_folder,glue::glue("Sample_correlation ",dataset_name,".pdf")))
    #plot_volcano(dep, contrast = "T0_vs_T24", label_size = 2, add_names = TRUE)
    # plot_heatmap(dep, type = "centered", kmeans = TRUE, 
    #              k = 6, col_limit = 4, show_row_names = FALSE,
    #              indicate = c("condition", "replicate"))
    # plot_single(dep, proteins = c("P05386","Q9H9B4"))
    # 
    # data_imp@assays@data@listData[[1]] %>% 
    #     as.data.frame() %>% 
    #     rownames_to_column("ProteinGroup") %>% 
    #     mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
    #     left_join(Interesting_proteins) %>% 
    #     subset(!is.na(Behaviour)) %>% 
    #     pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     mutate(Condition= factor(Condition)) %>% 
    #     group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
    #     mutate(DMSO = mean(c_across(contains("dmso")), na.rm = T)) %>% 
    #     mutate(across(where(is.numeric), ~.x-DMSO)) %>%
    #     dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
    #     ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID, alpha= Significant))+
    #     geom_line()+
    #     geom_point()+
    #     scale_alpha_manual(values = c(0.3,1))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == T))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == F))+
    #     theme(legend.position = "none") +
    #     facet_wrap("Behaviour")+ 
    #     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    #     ggtitle("Interesting DDR proteins Detected",
    #             "Significant - opaque, non-significant Transparent")
    # ggsave(here::here(output_folder,glue::glue("Known_Behaviour ",dataset_name,".pdf")), width = 20, height = 20)
    # data_norm@assays@data@listData[[1]] %>% 
    #     as.data.frame() %>% 
    #     rownames_to_column("ProteinGroup") %>% 
    #     mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
    #     pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     mutate(Condition= factor(Condition)) %>% 
    #     group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
    #     mutate(Samples_median = median(c_across(where(is.numeric)), na.rm = T)) %>%  
    #     mutate(across(where(is.numeric), ~.x-Samples_median)) %>%
    #     dplyr::select(-Samples_median) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
    #     subset(Significant == T) %>%
    #     ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID))+
    #     #geom_line()+
    #     geom_point()+
    #     #scale_alpha_manual(values = c(0.3,1))+
    #     #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == T))+
    #     #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == F))+
    #     theme(legend.position = "none") +
    #     facet_wrap("ID")+ 
    #     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    #     ggtitle(glue::glue("Significant - Proteins",dataset_name))
    # ggsave(here::here(output_folder,glue::glue("Significant_proteins ",dataset_name,".pdf")), width = 20, height = 20)
    # Significant_protein_list <- data.frame(ProteinGroup = dep@elementMetadata$name,
    #                                        dmso_X0 = dep@elementMetadata$dmso_vs_x0h_significant,
    #                                        T0_X24 = dep@elementMetadata$x0h_vs_x24h_significant,
    #                                        dmso_X240 = dep@elementMetadata$dmso_vs_x24h_significant) %>% 
    #     mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot)))
    # 
    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        #mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        #left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #left_join(Interesting_proteins) %>% 
        #subset(!is.na(Behaviour)) %>% 
        #pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
        #group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        #mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        #dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>% dplyr::select(-Significant) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #   mutate(Condition = str_remove_all(Condition,"_.")) %>% 
        # ungroup() %>% 
        #   group_by(ProteinGroup,Condition) %>% 
        # dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = T)) %>% 
        pivot_wider(names_from = "Condition",values_from = "Abundance") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
        ungroup %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        subset(!is.na(ID)) %>% 
        column_to_rownames("ID") %>%
        dplyr::select(where(is.numeric)) %>%
        mutate(Rowmean = rowMeans(.),
               across(where(is.numeric),~.x- Rowmean)) %>% 
        dplyr::select(-Rowmean) %>% 
        as.matrix() 
    paletteLength <- 20
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins, na.rm = T)/paletteLength, max(Significant_proteins, na.rm = T), length.out=floor(paletteLength/2)))
    
    png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
    pheatmap::pheatmap(Significant_proteins[,experimental_design_DIA$label %>% sort()],cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", cluster_rows = T,
                       # scale = "row",
                       main = glue::glue(dataset_name, " \nSignificant Proteins Normalised - imputted"),color=myColor, breaks=myBreaks)
    dev.off()
    # if(names(dev.cur()) != "null device"){dev.off()}
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
    Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
           # i = (dep@elementMetadata %>% names() %>% str_subset("diff"))[2]
        contrast <- str_remove_all(i,"_diff")
        print(contrast)
        conditions <- contrast %>% str_match("([:graph:]*)_vs_([:graph:]*)$") %>% .[2:3]
        non_missing_in_all_comparison <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames(),
                                           data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames()) %>% unique()
        Imputted <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames(),
                      data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames()) %>% unique()
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Uniprot = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
            subset(Uniprot %in% non_missing_in_all_comparison) %>% 
            
            mutate(Imputted_comparison = Uniprot %in% Imputted,
                   Single_Uniprot = Uniprot %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
            # left_join(Metabolic_proteins, by  = c("Single_Uniprot" = "Uniprot") ) %>% 
            # dplyr::rename(Metabolic_library = Behaviour) %>% 
            # left_join(Interesting_proteins, by = c("Single_Uniprot" = "Uniprot")) %>% 
            left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type), by  = c("Single_Uniprot" = "Uniprot")) #%>% 
            # mutate(Metabolic_library = if_else(is.na(Metabolic_library),"Non-Metabolic", Metabolic_library))
        volcano_df %>% ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = log2_FC, alpha = significant))+
            geom_point()+
            geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
            ggrepel::geom_label_repel(data = . %>% subset(significant == T))+
            annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
            
            ggtitle(glue::glue("Diff Present", contrast),
                    subtitle = dataset_name)
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_significant",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
        
        Comparisons_list[[i]] <- volcano_df
        

        # ego2 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T) ,
        #               OrgDb        = org.Hs.eg.db,
        #               ont          = "ALL",
        #               keyType = "UNIPROT",
        #               #nPerm        = 1000,
        #               minGSSize    = 50,
        #               maxGSSize    = 500,
        #               pvalueCutoff = 0.05,
        #               verbose      = FALSE)
        # if((ego2@result %>% nrow)>0){
        #     enrichplot::ridgeplot(ego2, showCategory = 68)+
        #         ggtitle(glue::glue(dataset_name," ALL-GSEA",i))
        #     ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," ALL-GSEA.pdf")), height = 20, width  = 15)}
        # # ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T),
        #               OrgDb        = org.Hs.eg.db,
        #               ont          = "BP",
        #               keyType = "UNIPROT",
        #               #nPerm        = 1000,
        #               minGSSize    = 50,
        #               maxGSSize    = 500,
        #               pvalueCutoff = 0.01,
        #               verbose      = FALSE)
        # if((ego3@result %>% nrow)>0){
        # 
        #     enrichplot::ridgeplot(ego3 , showCategory = 68)+
        #         ggtitle(glue::glue(dataset_name,"BP-GSEA",i))
        #     ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," BP-GSEA.pdf")), height = 20, width  = 15)}

        # gene_list <- volcano_df %>%
        #   dplyr::select(Single_Uniprot,log2_FC) %>%
        #   left_join(Human_hsa, by = c("Single_Uniprot" ="Uniprot")) %>%
        #   na.omit()%>%arrange(-log2_FC) %>%   pull(log2_FC,ID)
        # 
        # kk2 <- gseKEGG(geneList     = gene_list,
        #                organism     = 'hsa',
        #                minGSSize    = 10,
        #                pvalueCutoff = 0.05,
        #                verbose      = FALSE)
        # if((kk2@result %>% nrow)>0){
        # 
        #   enrichplot::ridgeplot(kk2, showCategory = 100)+
        #     ggtitle(glue::glue(dataset_name," KEGG ",i))
        #   ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," KEGG.pdf")), height = 20, width  = 15)}
        # 
        # 
        # mkk2 <- gseMKEGG(geneList = gene_list,
        #                  organism = 'hsa',
        #                  minGSSize = 10,
        #                  pvalueCutoff = 0.05)
        # if((mkk2@result %>% nrow)>0){
        # 
        #   enrichplot::ridgeplot(mkk2, showCategory = 68)+
        #     ggtitle(glue::glue(dataset_name," MKEGG ",i))
        #   ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," MKEGG.pdf")), height = 20, width  = 15)}
        # 
        # 
    }
    purrr::imap(.x = Comparisons_list, ~.x  %>% select(Uniprot, log2_FC, significant,ID) %>% rename(
        "{.y}_FC" := log2_FC,
        "{.y}_significant" := significant)) %>% 
        purrr::reduce(full_join) %>% 
        ggplot(aes(x = c_c_vs_c_uc_diff_FC,
                   y = n_c_vs_n_uc_diff_FC, label = ID))+geom_point()+
        ggrepel::geom_label_repel(data = . %>% subset((c_c_vs_c_uc_diff_FC>1.5 &n_c_vs_n_uc_diff_FC<(-1.5)) |
                                                          (c_c_vs_c_uc_diff_FC<(-1.5) &n_c_vs_n_uc_diff_FC>1.5)))+
        ggtitle("Nuclear and Cytoplasmic changes upon Confinement")
    ggsave(here::here(output_folder, glue::glue(i," ",dataset_name,"Combined_Cyto_Nucle_FC.pdf")), height = 20, width  = 15)

    data_imp@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle = dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_imputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
    data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle =  dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_unimputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
    # replicates <- experimental_design_DIA %>% group_by(condition) %>% sample_n(2) %>% ungroup() %>% pull(label,condition)
    # input_df <- data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot")
    # named_vector = replicates[1:2]
    # rep_1 = named_vector[1]
    # rep_2 = named_vector[2]
    # MA_replicates(rep_1,rep_2,input_df)
    # 
    # MA_replicates <-  function(named_vector,named_vector_2,input_df){
    #     print(named_vector)
    #     print(named_vector_2)
    #     
    #     input_df %>% mutate(
    #        "{unique(names(named_vector))}" :=  diff({{named_vector}},{{named_vector_2}}))
    # }
    # 
    # data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
    #     na.omit() %>% 
    #     mutate()
    #     #glimpse() %>% 
    #     pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
    #     mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
    #     group_by(Condition,Uniprot) %>% 
    #     summarise(Mean_Abundance = mean(Abundance, na.rm = T),
    #               Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
    #     ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
    #     geom_point()+
    #     facet_wrap("Condition")+
    #     ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
    #             subtitle =  dataset_name)
    # ggsave(here::here(output_folder,glue::glue("CV_conditions_imputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
    # 
    data_matrices <- list(Imputted = data_imp@assays@data@listData[[1]],
                          Unimputted = data_norm@assays@data@listData[[1]], 
                          DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff")))
    write.xlsx(data_matrices[1:2], here::here("Datasets","Processed",glue::glue(dataset_name, "matrices.xlsx")), overwrite = T,rowNames = T)
    write.xlsx(data_matrices[[3]], here::here("Datasets","Processed",glue::glue(dataset_name, "volcano_DFs.xlsx")), overwrite = T)
    
    return(data_matrices)
}



correp <- function(data_matrix,subset_selection){
    pacman::p_load(e1071,cluster,CORREP)
    #example from  https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/CORREP/inst/doc/CORREP.pdf
    #The format is = each column is a condition, each row is a Gene Replicate in that condiition and and all the replicates together
    
    # subset_selection = Significant_proteins
    # data_matrix <-Data_matrices$Imputted  
    data_matrix <- data_matrix[rownames(data_matrix) %in%  subset_selection, ]
    conditions = Data_matrices$Imputted  %>% colnames() %>%
        str_remove_all("_[:digit:]$") %>% unique()
    replicates <- Data_matrices$Imputted  %>% colnames() %>%
        str_match("_([:digit:]*)$") %>% .[,2] %>% as.numeric() %>% max()
    Pivotted_conditions <- map(.x= conditions, ~data_matrix %>% 
                                   dplyr::select(contains(.x)) %>% 
                                   rownames_to_column("Uniprot") %>% 
                                   
                                   pivot_longer(-Uniprot,names_to = "Replicate",values_to = .x) %>% 
                                   mutate(Replicate = str_remove_all(Replicate,paste0(.x,"_"))) %>% 
                                   unite(Uniprot,c(Uniprot,Replicate), sep = ", ")) %>% 
        purrr::reduce(left_join)
    
    d0.std <- apply(Pivotted_conditions %>% column_to_rownames("Uniprot"), 1, function(x) x/sd(x))
    input <- t(d0.std)
    M <- cor.balance(input, m=replicates, G=nrow(data_matrix)) 
    colnames(M) <- rownames(data_matrix)
    rownames(M) <- rownames(data_matrix)
    M_dist <- 1-M
    d <- as.dist(M_dist)
    g<- diana(d)
    clusters <-  NbClust::NbClust(method = "complete", distance = "euclidean", data = M, index =  "silhouette")$Best.partition %>% 
        enframe(name = "gene", value = "cluster")
    fuzziness <-  advclust::fuzzy.CM(M,K=clusters$cluster %>% max(),m=2,max.iteration=100,threshold=1e-5,RandomNumber=1234)
    Confident_items <-  fuzziness@member[(((fuzziness@member>0.5) %>% matrixStats::rowSums2()) ==1) & (((fuzziness@member>0.8) %>% matrixStats::rowSums2()) ==1),] %>% 
        rownames()
    clusters <- fuzziness@hard.label[Confident_items]
    M_confi <-  M[Confident_items,Confident_items]
    row_ha = ComplexHeatmap::rowAnnotation(clusters = as.character(clusters))
    
    
    list(correp_matrix = M,
         correp_matrix_conf = M_confi ,
         clusters = clusters,
         heatmap = ComplexHeatmap::Heatmap(M_confi, left_annotation = row_ha))
    
}

is_function = function (expr) {
    if (! is_assign(expr))
        return(FALSE)
    value = expr[[3]]
    is.call(value) && as.character(value[[1]]) == 'function'
}

function_name = function (expr){
    as.character(expr[[2]])}

is_assign = function (expr){
    is.call(expr) && as.character(expr[[1]]) %in% c('=', '<-', 'assign')}


produce_ipath_map <- function(df,
                              column_name = NULL){
    ### takes as input a dataframe produced which converts ratios to colours
    ### and the columns to map and produces iPath3 images 
    #df = For_ipath_comparisons
    #column_name = columns_for_mapping[1]
    
    Filtered_proteins <- df %>% 
        dplyr::select(all_of(c("V1","Width",column_name))) %>% 
        na.omit()
    
    #Filtered_proteins$Width <- paste("W",as.character(ntile(Filtered_proteins[,i], 30)),sep="")
    Filtered_proteins$Ipath <- paste(Filtered_proteins$V1, 
                                     Filtered_proteins %>% pull(column_name), 
                                     Filtered_proteins$Width)
    
    file_name <- here::here("Project_Output",
                      paste0(column_name,
                             "_ipath_all_enzymes_bins.tsv"))
    
    selections <- paste('selection=',
                        paste(paste(Filtered_proteins$Ipath,
                                    "%0A ",sep = ""),
                              collapse = " "),
                        collapse = " ")
    
    export_type <- "export_type=svg"
    
    file_ipath <- paste("./../Project_Output/",
                        column_name,"_ipath_all_enzymes_bins.svg",
                        sep="")
    
    ipath_command <- paste(paste("curl -d ",
                                 selections, " -d ", 
                                 export_type, 
                                 " https://pathways.embl.de/mapping.cgi -o ", 
                                 sep = '\"'),
                           file_ipath,
                           sep="")
    system(ipath_command)
    print(column_name)
    write_tsv(Filtered_proteins,file_name,col_names = FALSE)
    #pacman::p_load("rsvg")
    bitmap <- rsvg::rsvg_raw( here::here("Project_Output",
                                   paste0(column_name,"_ipath_all_enzymes_bins.svg")),
                              width = 3600)
    
    cowplot::ggdraw()+cowplot::draw_image(bitmap)+ cowplot::draw_figure_label(column_name,"top", size = 65)
    
    
}

calculateCoveredProtein_sav <- function (sample_id,proteinIDs, markerproteins) {
    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", 
                      "N4", "C1", "C2", "C3", "C4", "C5", "M1", "M2")
    color.code <- c("gold", "orange", "salmon", "tomato2", "grey90", 
                    "grey70", "grey50", "grey30", "lightblue", "aquamarine", 
                    "cyan", "deepskyblue2", "turquoise3", "burlywood4", 
                    "tan4")
    compartment.size <- c(358, 351, 252, 174, 192, 121, 231, 
                          198, 242, 132, 220, 215, 341, 69, 269)
    covered.proteins <- intersect(proteinIDs, markerproteins)
    if (length(covered.proteins) < 1) 
        warning("There is no overlap between marker proteins and data!")
    c.marker.df <- SubCellBarCode::markerProteins[covered.proteins, 
    ]
    coverageCompWise <- lapply(seq_len(length(compartments)), 
                               function(x) {
                                   temp.df <- c.marker.df[c.marker.df$Compartments == 
                                                              compartments[x], ]
                                   values <- list(Compartments = compartments[x], ColorCode = color.code[x], 
                                                  ProteinCoverage = 100 * ((dim(temp.df)[1])/compartment.size[x]))
                               })
    coverage.df <- as.data.frame(do.call("rbind", coverageCompWise))
    non.enriched.loc <- coverage.df[coverage.df$ProteinCoverage < 
                                        20, ]
    if (nrow(non.enriched.loc) == 1) {
        warning("There is not enough enrichment at: ", as.character(non.enriched.loc$Compartments), 
                "\nWe recommend you to perform the fractionation, again.")
    }
    else if (nrow(non.enriched.loc) > 1) {
        comp <- paste(as.character(non.enriched.loc$Compartments), 
                      collapse = ",")
        warning("There are not enough enrichments at: ", comp, 
                "\nWe recommend you to perform the fractionation!")
    }
    coverage.df$ProteinCoverage <- as.numeric(coverage.df$ProteinCoverage)
    coverage.df$Compartments <- as.character(coverage.df$Compartments)
    coverage.df$ColorCode <- as.character(coverage.df$ColorCode)
    ggplot(data = coverage.df, aes(x = coverage.df$Compartments, 
                                         y = coverage.df$ProteinCoverage)) + geom_bar(stat = "identity", 
                                                                                      fill = coverage.df$ColorCode) + 
        scale_x_discrete(limits = c(compartments)) + 
              theme_bw() + theme(text = element_text(size = 5), plot.title = element_text(hjust = 0.5), 
                                 axis.text.x = element_text(face = "bold", color = "black"), 
                                 axis.text.y = element_text(face = "bold", color = "black")) +
        ylim(0, 100)+
        
              labs(title = paste("Marker Protein Coverage Sample",sample_id), 
                   y = "% Protein Coverage", x = "Compartment")

}
