if (!require("pacman")) install.packages("pacman", repos = 'http://cran.us.r-project.org')
library(pacman)
#options(repos = getOption("repos")["CRAN"])
####MASCOT
#https://miktex.org/download
#if (!require("packrat")) install.packages("packrat", repos = 'http://cran.us.r-project.org')
#library(packrat)

### installing other packages
#pacman::p_load_gh("krlmlr/here" ,  dependencies = TRUE)
pacman::p_load("readr","tidyverse", "openxlsx", "stringr", "ontologyIndex", "biomaRt",
               "data.table","BiocManager","openxlsx","tidyverse","factoextra", "devtools",
               "proteusSILAC", "DescTools", "UpSetR", "VennDiagram")

if (!require("hpar")){
    BiocManager::install("hpar")
    library("hpar")
}
#needed for readGAF, and GoAnnotations, installable through BiocManager
# BiocManagerlibraries <- c("biomaRt","org.Hs.eg.db") #"clusterProfiler","pathview",
# 
# for (i in BiocManagerlibraries){
#     if (!requireNamespace(i, quietly = TRUE))
#         BiocManager::install(i)
#     library(i, character.only = T )

Process_Mascot <-  function(y){
    x <- read.xlsx(y)
    x <- x[,str_detect(names(x),"^Accession|Description|Symbol|Abundance:")]
    x[,4] <- as.numeric(x[,4])
    x[,5] <- as.numeric(x[,5])
    x$L_H_Ratio <- x[,4]/x[,5]
    x$Log2_LH <- log2(x$L_H_Ratio)
    x$Experiment <- str_match(y,"_(00.)-")[,2]
    x
    }

filenames <- list.files("./Project_Datasets/2020LD002", pattern="*.xlsx", full.names=TRUE)
ldf <- lapply(filenames, Process_Mascot)# %>% lapply()

# Common_SS_proteins <- Reduce(function(...) intersect(...), lapply(ldf[2:5], na.omit(as.vector("[")), 1)) %>%
#     unlist()

Heavy_SiLAC_present <- function(x){
    x$Accession[!is.na(x[,5])] %>% unlist()
}
Light_SiLAC_present <- function(x){
    x$Accession[!is.na(x[,4])] %>% unlist()
}

ldf2 <- Reduce(function(...) inner_join(..., by = "Accession"), lapply(ldf[2:5], commonly_detected_SS))



SuperSilac_comparison <- KendallW(ldf2[,-1],  test = T)

Heavy_SiLAC_presence_list <- lapply(ldf[2:5], Heavy_SiLAC_present)
Light_SiLAC_presence_list <- lapply(ldf, Light_SiLAC_present)

Common_SS_proteins <- Reduce(function(...) intersect(...), Heavy_SiLAC_presence_list)
commonly_detected_SS <- function(x){
    #x <- ldf[[3]]
    x <- x[!is.na(x[,5]) & x[,1] %in% Common_SS_proteins,c(1,5)]
    x[,2] <- rank(x[,2])
    x
    
}

names(Heavy_SiLAC_presence_list) <- paste("Sample", 2:5)
names(Light_SiLAC_presence_list) <- paste("Sample", 1:5)
upset(fromList(Heavy_SiLAC_presence_list), order.by = "freq")
dev.print(pdf, './Project_Output/SuperSilac_Comparison_Presence.pdf')

upset(fromList(Light_SiLAC_presence_list), order.by = "freq")
dev.print(pdf, './Project_Output/Light_Comparison_Presence.pdf')

SS_chrom_identified <- setdiff(Heavy_SiLAC_presence_list[[2]],Heavy_SiLAC_presence_list[[1]])
SS_chrom_identified_abundance <- summary(ldf[[3]]$`Abundance:.F3:.Heavy,.Control`[ldf[[3]]$Accession %in% SS_chrom_identified])
Other_SS_sample3_abundance <- summary(ldf[[3]]$`Abundance:.F3:.Heavy,.Control`[!(ldf[[3]]$Accession %in% SS_chrom_identified)])

SS_chrom_identified_abundance_Log <- summary(ldf[[3]]$Log2_LH[ldf[[3]]$Accession %in% SS_chrom_identified])
other_chrom_identified_abundance_Log <- summary(ldf[[3]]$Log2_LH[!(ldf[[3]]$Accession %in% SS_chrom_identified)])

Combining_dfs <- function(x){
    #x <- ldf[[2]]
    x <- x[,c(1,4,5,8)]
    colnames(x) <- c("Accession", "Light", "Heavy","Experiment")
    pivot_longer(x, c(Light, Heavy), names_to="Label", values_to = "Abundance" )
}

pivotted <- lapply(ldf,Combining_dfs) %>% bind_rows() %>% unique()
cor(pivotted$Experimentpivotted$Experiment == "Suit2" & pivotted$Label == "Light",)
# colnames(Maria_5M) <- c("MCF7_cytoplasm","MCF7_chromatome","NormalBreast_cytoplasm",'NormalBreast_chromatome')
# row.names(Maria_5M) <- openxlsx::read.xlsx("./Project_Datasets/Maria_5M_MS.xlsx")[,3]

# 
# Maria_5M[is.na(Maria_5M)] <- runif(length(Maria_5M[is.na(Maria_5M)]), min= Quantiles[1], max=Quantiles[2])
# Maria_5M <- as.matrix(Maria_5M)
# 
# res.dist <- dist(scale(log(t(Maria_5M))), method = "euclidean")
# res.hc <- hclust(d = res.dist, method = "ward.D2")
# fviz_dend(res.hc, cex = 0.5)

##### Checking cancer

#breast_cancer_set <- c("CCNB2", "FBXO5", 'KIF4A', 'MCM10', 'TPX2')


#pivotted <- pivot_longer(uninferred_Maria_5M, cols = colnames(uninferred_Maria_5M)[-1], names_to = "Sample", values_to="Abundance" )
Quantiles <- quantile(pivotted$Abundance,prob = seq(0, 1, length = 11), type = 5,na.rm = T)
HUMAN_9606 <- read_tsv("./Project_Datasets/HUMAN_9606_idmapping.dat",
                 col_names = FALSE) 

HUMAN_9606_idmapping <- HUMAN_9606 %>% filter(X2 == 'Gene_Name')

pivotted <- right_join(HUMAN_9606_idmapping[,-2],pivotted, by = c("X1" = "Accession"))
colnames(pivotted) <- c("Uniprot", "Gene", "Experiment","Label", "Abundance")


pivotted$Experiment <- dplyr::case_when(
    pivotted$Experiment == "001" ~ "Suit2",
    pivotted$Experiment == "002" ~ "HCT116",
    pivotted$Experiment == "003" ~ "MCF7",
    pivotted$Experiment == "004" ~ "MDAMB231",
    pivotted$Experiment == "005" ~ "HeLa",
    TRUE ~ as.character(pivotted$Experiment)
)
breast_cancer_set <- c("O95239", "Q9ULW0","P62805","P10109","P18206", "Q92783", "P84077", "P84243")
p <- subset(pivotted, Uniprot %in% breast_cancer_set) %>% 
    subset(!duplicated(.[,-2]))%>%
    ggplot( aes(x= Experiment, y = log(Abundance), colour = Label))+
    geom_col(aes(fill = factor(Gene)),position = position_dodge(width = 0.9), size =1 ) #aes(shape = Uniprot))
p + geom_hline(yintercept = log(Quantiles[1]), linetype="dashed", color = "black", size=0.5)+
    ggtitle("Subcompartment Markers for SuperSILAC and \n Chromatome (MCF7, MDAMB231, HeLa)")
ggsave("./Project_Output/SuperSilac_Markers_Confirmation_Proteins_BarChart.jpg")

POI <- c("P00338", "P12268","P43490","P13995","P07195", "P11586", "P04797", "O60885","P04406")
p <- subset(pivotted, Uniprot %in% POI) %>% 
    subset(!duplicated(.[,-2]))%>%
    ggplot( aes(x= Experiment, y = log(Abundance), colour = Label))+
    geom_col(aes(fill = factor(Gene)),position = position_dodge(width = 0.9))  
    #ylim(log(Quantiles[1]),log(max(pivotted$Abundance,na.rm =  T))) #aes(shape = Uniprot))
p + geom_hline(yintercept = log(Quantiles[1]), linetype="dashed", color = "black", size=0.5) +
    ggtitle("POI for SuperSILAC and \n Chromatome (MCF7, MDAMB231, HeLa)") + coord_cartesian(ylim = c(log(Quantiles[1]), log(Quantiles[10])+5))
ggsave("./Project_Output/SuperSilac_POI_BarChart.jpg")


x <- pivotted$Uniprot[!is.na(pivotted$Abundance)]
y <- subset(HUMAN_9606, X1 %in% x & X2 == "Ensembl")[,3]
HPA <- getHpa(unlist(y), hpadata = "hpaSubcellularLoc")  



for(j in c("Light","Heavy")){
    for(i in unique(pivotted$Experiment)){
        print(i)
        print(j)
        Uniprot_ids <- pivotted$Uniprot[!is.na(pivotted$Abundance) & pivotted$Label == j & pivotted$Experiment ==i]
        Ensembl_ids <- subset(HUMAN_9606, X1 %in% Uniprot_ids & X2 == "Ensembl")[,3]
        HPA2 <- getHpa(unlist(Ensembl_ids), hpadata = "hpaSubcellularLoc")$GO.id %>% str_split(";")%>% unlist()%>% table()%>%as.data.frame()
        HPA2 <- separate(HPA2,., c("Location", "GO"),sep = " \\(")
        ggplot(HPA2,aes(x = Location, y = Freq))+
            geom_col()+
            theme(axis.text.x  = element_text(angle=45, hjust = 1))+ 
            ggtitle(paste0("HPA GO terms ",i," ",j))
        ggsave(filename = paste0("./Project_Output/","HPA GO terms ",i," ",j,".jpg"))
    }
}
Nuclear_proteins <- HPA[str_detect(HPA$Enhanced, "Nucl") & !(str_detect(HPA$Enhanced, ";")),c(2,4)]

Nuclear_enrichment <- inner_join(Nuclear_proteins,pivotted, by = c("Gene.name" = "Gene")) %>% 
    pivot_wider(names_from = "Label", values_from = "Abundance") 

for(i in unique(Nuclear_enrichment$Experiment)){
    Nuclear_enrichment[Nuclear_enrichment$Experiment == i,]%>% 
        ggplot(aes(x = log10(Light), y = log10(Heavy), colour = Enhanced)) +
        geom_point() + 
        ggtitle(paste0("HPA Nucl Enhanced",i))
    ggsave(filename = paste0("./Project_Output/HPA Nucl Enhanced",i,".jpg"))
}


library(SubCellBarCode)

for(j in 3:4){
    for(i in unique(pivotted$Experiment)){
        Heavy_or_light <- if_else(j == 3, "light", "heavy")
        x <-  pivotted[pivotted$Experiment == i,-2] %>% unique() %>% pivot_wider(names_from = "Label", values_from = "Abundance") %>%
            as.data.frame()
        #nam <- pivotted$Uniprot
        #x <- left_join(x, HUMAN_9606[HUMAN_9606$X2 == "GI",],by = c("Uniprot" = "X3")) %>% na.omit()
        #x <- x[!duplicated(x[ , "X1"]),]
        rownames(x) <- x$Uniprot
        #assign(nam,x)
        df <- convert2symbol(df = x, id = "UNIPROT")
        c.prots <- calculateCoveredProtein(proteinIDs = rownames(subset(df, !is.na(df[,j]))),markerproteins = markerProteins[,1])
        ggsave(filename = paste("./Project_Output/Fractionation_Markers",i,j,".jpg", sep = "_"))
    }
}
chromatome_genes <- markerProteins$Proteins[markerProteins$Compartments == "N3"]
library(ggpubr)
# Grouped Scatter plot with marginal density plots
pivotted_wider <- pivotted %>% pivot_wider(names_from = "Label", values_from = "Abundance") %>% .[!(is.na(.$Light) & is.na(.$Heavy)),]
pivotted_wider_cor <- pivotted_wider %>% na.omit()

scatter_plot <- subset(pivotted_wider_cor, Experiment == "Suit2") %>% ggplot(aes(x = log10(Light),  y= log10(Heavy)))
scatter_plot + geom_point() + labs(x = "Log10 - Light", y = "Log10 - Heavy") + geom_smooth(method="lm") + ggtitle("Suit2_light_vs_Heavy")
ggsave(filename = './Project_Output/Suit2_light_vs_Heavy_scatter.png')

cor(rank(log10(pivotted_wider_cor$Light[pivotted_wider_cor$Experiment == "Suit2"])), 
    rank(log10(pivotted_wider_cor$Heavy[pivotted_wider_cor$Experiment == "Suit2"])),
    "complete.obs", method = "spearman")
venn.diagram(
    x = list(pivotted_wider$Uniprot[pivotted_wider$Experiment == "Suit2" & !is.na(pivotted_wider$Light)],
             pivotted_wider$Uniprot[pivotted_wider$Experiment == "Suit2" & !is.na(pivotted_wider$Heavy)]),
    category.names = c("Suit2 Light" , "Suit2 Heavy" ),
    filename = './Project_Output/Suit2_light_vs_Heavy.png',
    output=TRUE
)


pivotted_wider$Light[is.na(pivotted_wider$Light)] <- quantile(pivotted_wider$Light, na.rm = T)[2]
pivotted_wider$Log10Light <- log10(pivotted_wider$Light)
pivotted_wider$Heavy[is.na(pivotted_wider$Heavy)] <- quantile(pivotted_wider$Heavy, na.rm = T)[2]
pivotted_wider$Log2Ratio <- log2(pivotted_wider$Light/pivotted_wider$Heavy)
pivotted_wider$Log10Ratio <- log10(pivotted_wider$Light/pivotted_wider$Heavy)
pivotted_wider$Compartment <- if_else(pivotted_wider$Gene %in% chromatome_genes,"Chromatin","Non-Chromatin")

for(i in unique(pivotted_wider$Experiment)){
    ggscatterhist(
        pivotted_wider[pivotted_wider$Experiment == i,], x = "Log2Ratio", y = "Log10Light",
        color = "Compartment", size = 1, alpha = 0.6,
        palette = c("#00AFBB", "#FC4E07"),
        margin.params = list(fill = "Compartment", color = "black", size = 0.2),
        title = i)
    ggsave(filename = paste0("./Project_Output/N3_Chrom_markers_Enrichment", i,".jpg"))
}

for(i in unique(pivotted_wider$Experiment)){
    pivotted_wider[pivotted_wider$Experiment == i,] %>%  ggplot(aes(x=Log2Ratio)) + 
        geom_histogram(color="black", fill="white") + 
        ggtitle(paste0("Log2 Ratio Light vs Heavy ", i))

    ggsave(filename = paste0("./Project_Output/Log2Ratio Histogram", i,".jpg"))
}

##uninferred Maria 5M from Maria_5M test
Last_MCF7_Chrom <- uninferred_Maria_5M[,c(1,3)] %>% na.omit()
venn.diagram(
    x = list(pivotted_wider$Uniprot[pivotted_wider$Experiment == "MCF7" & !is.na(pivotted_wider$Light)],
             Last_MCF7_Chrom$Uniprot),
    category.names = c("MCF7 \nCHAPS Chrom" , "Last \nMCF7 Chrom" ),
    filename = './Project_Output/MCF7_Chrom_New_OLD.png',
    output=TRUE
)
Cor_df <- inner_join(Last_MCF7_Chrom,pivotted_wider[pivotted_wider$Experiment == "MCF7" & !is.na(pivotted_wider$Light),]) 
cor(rank(as.numeric(Cor_df$MCF7_chromatome)), rank(Cor_df$Light), method = "spearman")
# Euchrom <- c("Q03164", "Q9UKN8","Q6ZRS2","Q8WTS6","Q9UMN6", "O95251", "Q92794", "Q92793","Q86X55",
#              "P21675", "P15336", "Q99873", 'O14929')
# p <- subset(pivotted, Uniprot %in% Euchrom) %>% 
#     subset(!duplicated(.[,-2]))%>%
#     ggplot( aes(x= Sample, y = log(Abundance), colour = Gene))+
#     geom_col(aes(fill = factor(Gene)),position = position_dodge(width = 0.9)) #aes(shape = Uniprot))
# p + geom_hline(yintercept = log(Quantiles[1]), linetype="dashed", color = "black", size=0.5)
# ggsave("./Project_Output/Maria5M_Sample_EUchrom_BarChart.jpg")
# 
# Heterochrom <- c("Q9NQR1", "Q96KQ7","Q15910","Q9H9B1","Q9H5I1", "P26358", "Q15022","O75530",
#                  "Q9Y232","06587","O43463", "Q4FZB7", "Q9Y6K1","Q86Y97")
# p <- subset(pivotted, Uniprot %in% Heterochrom) %>% 
#     subset(!duplicated(.[,-2]))%>%
#     ggplot( aes(x= Sample, y = log(Abundance), colour = Gene))+
#     geom_col(aes(fill = factor(Gene)),position = position_dodge(width = 0.9)) #aes(shape = Uniprot))
# p + geom_hline(yintercept = log(Quantiles[1]), linetype="dashed", color = "black", size=0.5)
# ggsave("./Project_Output/Maria5M_Sample_Heterochrom_BarChart.jpg")

### http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/proteus.html#clustering
###MAXQUANT

#install.packages("proteus")

#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("limma")
#install.packages("devtools")
#devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
library("proteus")

devtools::install_github("bartongroup/proteusSILAC")
library(proteusSILAC)
data(proteusSILAC)


measCols <- list(
    HL = "Ratio H/L"
)

eviCols <- list(
    sequence = 'Sequence',
    modified_sequence = 'Modified sequence',
    modifications = 'Modifications',
    protein_group = 'Proteins',
    protein = 'Leading razor protein',
    experiment = 'Experiment',
    charge = 'Charge',
    reverse = 'Reverse',
    contaminant = 'Potential contaminant'
)


evidenceFile <- "./SILAC_02092020/combined/txt/evidence.txt"
#testing <- read.delim(evidenceFile)
evi <- readEvidenceFile(evidenceFile, measure.cols=measCols, data.cols=eviCols, zeroes.are.missing=FALSE)
metadataFile <- "SILAC_02092020/combined/txt/meta.txt"
meta <- read.delim(metadataFile, header=TRUE, sep="\t")
pepdat <- makePeptideTable(evi, meta, measure.cols=measCols, aggregate.fun=aggregateMedian, experiment.type="SILAC")
summary(pepdat)
plotCount(pepdat)
prodat <- makeProteinTable(pepdat, aggregate.fun=aggregateMedian)
summary(prodat)
prodat.norm <- normalizeData(prodat)
plotSampleDistributions(prodat, fill="replicate") + labs(title="Before")
plotSampleDistributions(prodat, fill="replicate") + labs(title="After")
# luni <- lapply(as.character(prodat$proteins), function(prot) {
#     if(grepl("sp\\|", prot)) {
#         uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
#         c(prot, uniprot)
#     }
# })
# ids <- as.data.frame(do.call(rbind, luni))
# names(ids) <- c("protein", "uniprot")
annotations <- fetchFromUniProt(prodat$proteins, verbose=TRUE)
annotations.id <- merge(ids, annotations, by.x="uniprot", by.y="id")
annotations.id <- unique(annotations.id)
prodat.med <- annotateProteins(prodat.med, annotations.id)


evi <- readEvidenceFile("./Project_Datasets/2020LD000_MAGU_20200617_evidence.txt")
##this depends on experiment
meta <- data.frame(experiment = 1:4, measure = "Intensity", 
                   sample = c("MCF7_Cyto", "MCF7_Chrom", "Norm_Cyto", "Norm_Chrom"),
                   condition = 1:4)
pepdat <- makePeptideTable(evi, meta)
plotCount(pepdat)
plotDistanceMatrix(pepdat) #peptidedistancematrix
ggsave("./Project_Output/Peptide_Distance_Matrix.jpg")
plotClustering(pepdat)
ggsave("./Project_Output/Peptide_Dendrogram.jpg")
plotPCA(pepdat)
ggsave("./Project_Output/Peptide_PCA.jpg")

prodat <- makeProteinTable(pepdat)
summary(prodat)

prodat.med <- normalizeData(prodat)
plotSampleDistributions(prodat, title="Not normalized", fill="condition", method="violin")
ggsave("./Project_Output/Protein_Sample_Distributions.jpg")
plotSampleDistributions(prodat.med, title="Median normalization", fill="condition", method="violin")
ggsave("./Project_Output/Norm_Protein_Sample_Distributions.jpg")
plotClustering(prodat.med)
ggsave("./Project_Output/Protein_Dendro_Maxquant.jpg")

luni <- lapply(as.character(prodat$proteins), function(prot) {
    if(grepl("sp\\|", prot)) {
        uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
        c(prot, uniprot)
    }
})
ids <- as.data.frame(do.call(rbind, luni))
names(ids) <- c("protein", "uniprot")
annotations <- fetchFromUniProt(ids$uniprot, verbose=TRUE)
annotations.id <- merge(ids, annotations, by.x="uniprot", by.y="id")
annotations.id <- unique(annotations.id)
prodat.med <- annotateProteins(prodat.med, annotations.id)

## SEE 5.3.1 Aggregator function We recommend using this function with SILAC data
 
#Read about limma in the website http://www.compbio.dundee.ac.uk/user/mgierlinski/proteus/proteus.html#clustering
#res <- limmaDE(prodat.med, sig.level=0.05) this is only for pairs, for more complicated, check http://bioconductor.org/packages/release/bioc/html/limma.html
#find which proteins are only detected in 1 conditions across all replicates
#only.A <- which(prodat$detect$`A` & !prodat$detect$B)
#only.B <- which(!prodat$detect$`A` & prodat$detect$B)

### ARTMS
#BiocManager::install("artMS")
library(artMS)
artmsWriteConfigYamlFile(config_file_name = "my_config.yaml")
artmsQualityControlEvidenceBasic(evidence_file = "./Project_Datasets/evidence.txt", 
                                 keys_file = "./Project_Datasets/keys.txt" )
artmsQuantification(yaml_config_file = "my_config.yaml")

###SubCellular
#BiocManager::install("SubCellBarCode")
library(SubCellBarCode)
df <- convert2symbol(df = Maria_5M, id = "UNIPROT")
#Maria_5M_Subcell <- Maria_5M %>% mutate()
Sample_names <- c("MCF7_cyto", "MCF7_chrom", "Norm_cyto", "Norm_chrom" )
for(i in 1:4){
    c.prots <- calculateCoveredProtein(proteinIDs = rownames(subset(df, !is.na(df[,i]))),markerproteins = markerProteins[,1])
    ggsave(filename = paste("./Project_Output/Fractionation_Markers_",Sample_names[i],".jpg", sep = ""))
}
Sample_old_names <- c("old_231", "old_BT474", "old_MCF7", "old_SKBR3","old_T47D" )
Laura1_Chrom <- read.xlsx("./Project_Datasets/MS protein ungrouped mean values.xlsx")[,c(4,33:37)]
rownames(Laura1_Chrom) <- Laura1_Chrom[,1]
df_2 <- convert2symbol(df = Laura1_Chrom[2:6], id = "UNIPROT")
for(i in 1:5){
    c.prots <- calculateCoveredProtein(proteinIDs = rownames(subset(df_2, !is.na(df_2[,i]))),markerproteins = markerProteins[,1])
    ggsave(filename = paste("./Project_Output/Fractionation_Markers_",Sample_old_names[i],".jpg", sep = ""))
}

MetaPhOrs <- read.csv("./Project_Datasets/MetaPhOrs.csv")[,1:2] %>%  
    mutate(Mus.musculus = str_trim(str_remove_all(.[,1],pattern = "Â"))) %>% 
    mutate(Homo.sapiens = str_trim(str_remove_all(.[,2],pattern = "Â"))) 
    
Edu_dataset <- read.xlsx("./Project_Datasets/Edu_Chromatome.xlsx")[-c(1:4),c(1,9:22)]
rownames(Edu_dataset) <- Edu_dataset[,1]
Edu_dataset <- data.matrix(Edu_dataset[,-1], rownames.force = NA) %>% 
    subset(rowSums(.[,(ncol(.)-8):ncol(.)],na.rm = T)>0 & rowSums(.[,1:7] == 0)>5)
Edu_dataset_human_homo <- subset(MetaPhOrs, Mus.musculus %in% rownames(Edu_dataset))[,2]
df_3 <- Edu_dataset[1:length(Edu_dataset_human_homo),] ###incorrectly assembled names
rownames(df_3) <- Edu_dataset_human_homo
df_3 <- convert2symbol(df = df_3, id = "UNIPROT")
c.prots <- calculateCoveredProtein(proteinIDs = rownames(df_3),markerproteins = markerProteins[,1])
ggsave(filename = paste("./Project_Output/Fractionation_Markers_","Edu_control_removed",".jpg", sep = ""))

Dutta <- read.xlsx("./Project_Datasets/Dutta_Chromatome_Hypoxia_Protein_Level.xlsx",sheet =2)[-c(1:13),-c(1:5)] %>% na.omit()
rownames(Dutta) <- Dutta[,1]
df_4 <- convert2symbol(df = Dutta, id = "UNIPROT")
c.prots <- calculateCoveredProtein(proteinIDs = rownames(df_4),markerproteins = markerProteins[,1])
ggsave(filename = paste("./Project_Output/Fractionation_Markers_","Dutta_Hypoxia_Chrom",".jpg", sep = ""))

Sdelci_WCE_trial <- read.xlsx("./Project_Datasets/Sdelci_Trial_data.xlsx",sheet =5)[-c(1:2),-c(1:2)]
rownames(Sdelci_WCE_trial) <- Sdelci_WCE_trial[,1]
df_5 <- convert2symbol(df = Sdelci_WCE_trial, id = "UNIPROT")
c.prots <- calculateCoveredProtein(proteinIDs = rownames(df_5),markerproteins = markerProteins[,1])
ggsave(filename = paste("./Project_Output/Fractionation_Markers_","Sdelci_WCE_nonCross",".jpg", sep = ""))

Sdelci_MCF7noCross_trial <- read.xlsx("./Project_Datasets/Sdelci_Trial_data.xlsx",sheet =3)[-c(1:2),-c(1:2)] %>% na.omit()
rownames(Sdelci_MCF7noCross_trial) <- Sdelci_MCF7noCross_trial[,1]
df_6 <- convert2symbol(df = Sdelci_MCF7noCross_trial, id = "UNIPROT")
c.prots <- calculateCoveredProtein(proteinIDs = rownames(df_6),markerproteins = markerProteins[,1])
ggsave(filename = paste("./Project_Output/Fractionation_Markers_","Sdelci_MCF7noCross_trial",".jpg", sep = ""))

Sdelci_MCF7Cross_trial <- read.xlsx("./Project_Datasets/Sdelci_Trial_data.xlsx",sheet =2)[-c(1:2),-c(1:2)] %>% na.omit()
rownames(Sdelci_MCF7Cross_trial) <- Sdelci_MCF7Cross_trial[,1]
df_7 <- convert2symbol(df = Sdelci_MCF7Cross_trial, id = "UNIPROT")
c.prots <- calculateCoveredProtein(proteinIDs = rownames(df_7),markerproteins = markerProteins[,1])
ggsave(filename = paste("./Project_Output/Fractionation_Markers_","Sdelci_MCF7Cross_trial",".jpg", sep = ""))

###Venn Diagrams
#install.packages("eulerr")
#install.packages("VennDiagram")
library("VennDiagram")
MCF_Old_New_chrom_Venn <- VennDiagram::get.venn.partitions(list(MCF7_new=rownames(subset(df, !is.na(df[,2]))), MCF7_old= rownames(subset(df_2, !is.na(df[,3]))) ))
grid.newpage()
VennDiagram::venn.diagram(list(MCF7_new=rownames(subset(df, !is.na(df[,2]))), MCF7_old= rownames(subset(df_2, !is.na(df[,3]))) ), 
                                          "./Project_Output/Old_vs_new_protocol_MCF7_Venn_chrom.jpg")
MCF_Old_New_chrom_Venn <- list(Venn = MCF_Old_New_chrom_Venn, 
                               Both_MCF7_new_MCF7_old = MCF_Old_New_chrom_Venn$..values..$`1`,
                               Only_old_protocol = MCF_Old_New_chrom_Venn$..values..$`2`,
                               Only_new_protocol = MCF_Old_New_chrom_Venn$..values..$`3`)
write.xlsx(MCF_Old_New_chrom_Venn,"./Project_Output/Old_vs_new_protocol_MCF7_Venn_chrom.xlsx")
#MaxQuant vs Mascot Venn

MaxQuant_With_Uniprot <- left_join(prodat.med[["stats"]], prodat.med[["annotation"]][,1:2], by = c("id" = "protein"))
VennDiagram::venn.diagram(list(MCF7_cyto_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$MCF7_cytoplasm),1],
                               MCF7_cyto_Maxquant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 1 & !is.na(MaxQuant_With_Uniprot$mean),6], 
                               MCF7_chrom_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$MCF7_chromatome),1],
                               MCF7_chrom_MaxQuant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 2 & !is.na(MaxQuant_With_Uniprot$mean),6]),
                             "./Project_Output/MCF7_MaxQuant_Mascot.tiff")

VennDiagram::venn.diagram(list(Norm_cyto_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$NormalBreast_cytoplasm),1],
                               Norm_cyto_Maxquant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 3 & !is.na(MaxQuant_With_Uniprot$mean),6], 
                               Norm_chrom_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$NormalBreast_chromatome),1],
                               Norm_chrom_Maxquant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 4 & !is.na(MaxQuant_With_Uniprot$mean),6]),
                          "./Project_Output/Norm_MaxQuant_Mascot.tiff")

VennDiagram::venn.diagram(list(MCF7_cyto_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$MCF7_cytoplasm),1],
                               Norm_cyto_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$NormalBreast_cytoplasm),1],
                               MCF7_chrom_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$MCF7_chromatome),1],
                               Norm_chrom_Mascot=uninferred_Maria_5M[!is.na(uninferred_Maria_5M$NormalBreast_chromatome),1]),
                          "./Project_Output/Mascot_Sample_Venn.tiff")

VennDiagram::venn.diagram(list(MCF7_cyto_MaxQuant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 1 & !is.na(MaxQuant_With_Uniprot$mean),6],
                               Norm_cyto_MaxQuant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 3 & !is.na(MaxQuant_With_Uniprot$mean),6], 
                              MCF7_chrom_MaxQuant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 2 & !is.na(MaxQuant_With_Uniprot$mean),6],
                              Norm_chrom_MaxQuant=MaxQuant_With_Uniprot[MaxQuant_With_Uniprot$condition == 4 & !is.na(MaxQuant_With_Uniprot$mean),6]),
                          "./Project_Output/MaxQuant_Sample_Venn.tiff")



####ipath3
Metabolic_ids <- as.vector(read.xlsx("./Project_Datasets/uniprot_hsa01100.xlsx")[,2])#Metabolic Uniprots
Metabolic_detected<- read_tsv("./Project_Datasets/FINAL_Valid_terms_for_ipath.tsv")


For_ipath_comparisons <- left_join(Metabolic_detected[,1:3],uninferred_Maria_5M, by = "Uniprot")

#Uniprots_alternative_ids_metabolic$color <- color.scale(Uniprots_alternative_ids_metabolic$Log2_Abundance,c(0,1,1),c(1,1,0),0)
For_ipath_comparisons$Width <- 100 #paste("W",as.integer(-log10(For_ipath_comparisons$P_value)*10),sep="")
#For_ipath_comparisons <- For_ipath_comparisons[order(For_ipath_comparisons$Log2_Abundance),]

#rbPal <- colorRampPalette(c('blue','grey','red'))
For_ipath_comparisons$Col <- "green"#rbPal(10)[as.numeric(cut(For_ipath_comparisons$Log2_Abundance,breaks = 10))]  

# for(i in levels(as.factor(For_ipath_comparisons$Comparisons))){
#     #i <- levels(as.factor(For_ipath_comparisons$Comparisons))[1]
#     For_ipath_unique_comp <- For_ipath_comparisons[For_ipath_comparisons$Comparisons == i,]
#     filename <- paste(i,"Ipath3_dataset.tsv")
#     write_tsv(For_ipath_unique_comp[For_ipath_unique_comp$P_value<0.05,c("V1","Width","Col")],filename,col_names = FALSE)
# }
# For_ipath_core_all_abundances <- subset(Metabolic_detected[,1:3],Uniprot %in% unlist(read_tsv("Common_core_all_abundances.tsv")))
# For_ipath_core_all_abundances$Col <- "#72903F"
# For_ipath_core_all_abundances$Width <- "W20"
# write_tsv(For_ipath_core_all_abundances[,c("V1","Width","Col")],"ipath_common_core_all_abundances.tsv",col_names = FALSE)
# 
# Common_core_stable<- read_tsv("Common_core_stable.tsv")
# For_ipath_core_stable <- na.omit(left_join(Metabolic_detected[,1:3],Common_core_stable, by = "Uniprot"),cols = 'Combined_Mean_Uninferred_abundance')
# For_ipath_core_stable$Col <- "#72903F"
# For_ipath_core_stable$Width <- paste("W",as.character(ntile(For_ipath_core_stable$Combined_Mean_Uninferred_abundance, 30)),sep="")
# For_ipath_core_stable <- For_ipath_core_stable[,c("V1","Width","Col")]
# write.table(For_ipath_core_stable,"For_ipath_core_stable.txt",row.names = FALSE,col.names = FALSE, quote = FALSE)
# 
# #not useful as most of them will have very high p-values
# unstable_core <- setdiff(For_ipath_core_all_abundances[,1],For_ipath_core_stable[,1])

# #ranking the proteins in 30 bins according to abundance
# all_proteins_uninferred_NA_removed_meaned <- read_tsv("all_proteins_uninferred_NA_removed_meaned.tsv")
# ipath_all_proteins_uninferred_NA_removed_meaned <- inner_join(Metabolic_detected[,1:3],all_proteins_uninferred_NA_removed_meaned,by = "Uniprot")
# for (i in unique(unlist(ipath_all_proteins_uninferred_NA_removed_meaned[,"Cell_line"]))){
#     i <- unique(unlist(ipath_all_proteins_uninferred_NA_removed_meaned[,"Cell_line"]))[1]
#     print(i)
#     Filtered_proteins <-unique(ipath_all_proteins_uninferred_NA_removed_meaned[ipath_all_proteins_uninferred_NA_removed_meaned$Cell_line == i,])
#     Filtered_proteins$Col <- "#72903F"
#     Filtered_proteins$Width <- paste("W",as.character(ntile(Filtered_proteins$Mean_Uninferred_abundance, 30)),sep="")
#     Filtered_proteins <- Filtered_proteins[,c("V1","Width","Col")]
#     file_name <- paste(i,"ipath_all_enzymes_bins_30.tsv",sep="_")
#     write_tsv(Filtered_proteins,file_name,col_names = FALSE)
#     
# }

# all_proteins_uninferred_NA_removed_meaned <- read_tsv("all_proteins_uninferred_NA_removed_meaned.tsv")
# ipath_all_proteins_uninferred_NA_removed_meaned <- inner_join(Metabolic_detected[,1:3],all_proteins_uninferred_NA_removed_meaned,by = "Uniprot")
# top_50 <- ipath_all_proteins_uninferred_NA_removed_meaned
#top_50 <- ipath_all_proteins_uninferred_NA_removed_meaned[ipath_all_proteins_uninferred_NA_removed_meaned$Mean_Uninferred_abundance>median(ipath_all_proteins_uninferred_NA_removed_meaned$Mean_Uninferred_abundance,na.rm = TRUE),]
# hist(ipath_all_proteins_uninferred_NA_removed_meaned$Mean_Uninferred_abundance)
# hist(top_50$Mean_Uninferred_abundance)
# top_50$Width <- paste("W",as.character(ntile(top_50$Mean_Uninferred_abundance, 30)),sep="")
file_ipath_names <- vector()
for (i in 4:7){
    #i <- unique(unlist(ipath_all_proteins_uninferred_NA_removed_meaned[,"Cell_line"]))[3]
    print(i)
    Filtered_proteins <- subset(For_ipath_comparisons,!is.na(For_ipath_comparisons[,i]))
    Filtered_proteins$Col <- "#72903F"
    Filtered_proteins$Width <- paste("W",as.character(ntile(Filtered_proteins[,i], 30)),sep="")
    Filtered_proteins$Ipath <- paste(Filtered_proteins$V1, Filtered_proteins$Col, Filtered_proteins$Width)
    #Filtered_proteins <- Filtered_proteins[,c("V1","Width","Col")]
    file_name <- paste("./Project_Output/",colnames(Filtered_proteins)[i],"_ipath_all_enzymes_bins.tsv",sep="")
    selections <- paste('selection=',paste(paste(Filtered_proteins$Ipath,"%0A ",sep = ""),collapse = " "),collapse = " ")
    export_type <- "export_type=svg"
    
    file_ipath <- paste("./Project_Output/",colnames(Filtered_proteins)[i],"_ipath_all_enzymes_bins.svg",sep="")
    file_ipath_names <- append(file_ipath_names,file_ipath)
    ipath_command <- paste(paste("curl -d ",selections, " -d ", export_type, " https://pathways.embl.de/mapping.cgi -o ", sep = '\"'),file_ipath,sep="")
    system(ipath_command)
    write_tsv(Filtered_proteins,file_name,col_names = FALSE)
    
}
write(file_ipath_names,'file_ipath_names.tsv')
