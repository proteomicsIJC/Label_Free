### Load libraries and data ----
# libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#source("log2Transf.R")
#source("proteinGroupsCleanner.R")
library(ggplot2)
library(ggfortify)
library(limma)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(janitor)
library(dplyr)
library(tidyr)
library(mice)
library(reshape2)
library(openxlsx)
# data
# maxquant_data
maxquant <- read.csv2("./raw_data/proteinGroups.txt", header = T, check.names = "F", sep = "\t", dec = ".")

# Meta data
meta_data <- read.csv2("./raw_data/meta_data.csv", header = T, check.names = T, sep = ",")

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

# Check if folders exist, if NO create them
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))
#----

### Function deffinition----
use_this_intensity <- function(dataset,pattern){
  intensities <- grep(pattern = pattern, colnames(dataset))
  major <- grep(pattern = "Majority protein IDs", colnames(dataset))
  dataset <- dataset[,c(c(major),
                        c(intensities))]
  
  colnames(dataset) <- gsub(x = colnames(dataset), 
                            pattern = pattern, 
                            replacement = "")
  
  colnames(dataset) <- gsub(x = colnames(dataset),
                            pattern = " ",
                            replacement = "")
  dataset <- clean_names(dataset)
  dataset$accession <- dataset$majorityprotein_i_ds
  dataset <- dataset %>%
    select(accession, everything())
  dataset <- subset(dataset, select = -c(majorityprotein_i_ds))
  dataset$protein_group <- rownames(dataset)
  dataset <- dataset %>% 
    select(protein_group, everything())
  return(dataset)
}

zero_to_NA <- function(patterns,dataset) {
  print("Zero values have been transformed to NAs")
  for (i in 1:length(patterns)){
    searching <- grep(pattern = patterns[[i]], names(dataset))
    dataset[searching][dataset[searching] == 0] <- NA}
  return(dataset)
}

removecontaminants <- function(dataset, contaminants, accession_name) {
  print(paste0("Removing contaminants from the data using contaminants file"))
  cn <- colnames(dataset)
  #cont <- readLines(contaminants)
  cont <- cont[grep(">", cont)]
  cont <- substr(cont, 2, 7)
  cont <- cont[-which(cont %in% c("_",".","ENSEMB","REFSEQ","Strept","H-INV:"))]
  cleaned_dataset <- dataset[-c(unlist(lapply(cont, function(x) grep(x, dataset[,grep(accession_name, colnames(dataset))])))), ]
  return(cleaned_dataset)
}

log2_to_pattern <- function(patterns,dataset){
  print("log 2 transformation of the intensity values")
  print("Note that, negative values after transformation have been assigned to NA")
  cat("log2 transformation",file = "./results/used_parameters.txt",sep = "\n",append = T)
  cat(paste0("Intensity values have been log2 transformed and negative values after transformation have been assigned to NA"), file = "./results/used_parameters.txt", sep = "\n", append = T)
  cat(paste0(rep("_",50), collapse = ""), file = "./results/used_parameters.txt",append = T, sep = "\n")
  for (i in 1:length(patterns)){
    searching <- grep(pattern = patterns[[i]], names(dataset))
    dataset[searching] <- log2(dataset[searching])
    dataset[searching][dataset[searching] <= 0] <- NA}
  return(dataset)
}

remove_samp <- function(dataset, samples=NULL){
  if(!is.null(samples)){
    print("Removing indicated samples")
    for (i in 1:length(samples)){
      dataset <- subset(dataset, sample_name!=samples[i])}
    cat("remove_samp", file = "./results/used_parameters.txt",append = T, sep = "\n")
    cat(paste0("The following samples have been removed due to low quality of the data: ",samples), file = "./results/used_parameters.txt",append = T, sep = "\n")
    cat(paste0(rep("_",50), collapse = ""), file = "./results/used_parameters.txt",append = T, sep = "\n")
  }
  if (is.null(samples)){
    print("No will be removed")
    dataset <- dataset
    cat("remove_samp",file = "./results/used_parameters.txt",append = T, sep = "\n")
    cat("No samples have been removed", file = "./results/used_parameters.txt",append = T, sep = "\n")
    cat(paste0(rep("_",50), collapse = ""), file = "./results/used_parameters.txt",append = T, sep = "\n")
  }
  return(dataset)
}

tim <- function(impute,dataset,NAs_prop){
  # PGs that have some NA
  for_nas <- dataset
  cols <- c("protein_group")
  for_nas <- for_nas %>%
    group_by(across(all_of(cols))) %>%
    mutate(count_NA = sum(is.na(normalized_intensity))) %>%
    ungroup()
  for_nas <- for_nas %>%
    filter(count_NA > 0)
  pg_with_nas <<- for_nas
  
  # Only NAs values
  only_nas <<- dataset[!complete.cases(dataset),]
  n <- length(unique(dataset$sample_name))
  if (impute == "no"){
    print("No missing values imputation is performed, all PGs with more than 0 NA values will be removed")
    cols <- c("protein_group")
    dataset <- dataset %>%
      group_by(across(all_of(cols))) %>%
      mutate(count_NA = sum(is.na(normalized_intensity))) %>%
      ungroup()
    
    # Filter pg with less than NAs_prop
    dataset <- dataset %>%
      filter(count_NA == 0)
    
    # Filter pg with less than NAs_prop
    dataset <- dataset %>% 
      group_by(across(all_of(cols))) %>%
      filter(n() >= length(unique(sample_name))) %>%
      ungroup()
    imputed_data <- dataset
    cat("tim - for NA imputation",file = "./results/used_parameters.txt",sep = "\n", append = T)
    cat(paste0("impute == no"), file = "./results/used_parameters.txt", sep = "\n", append = T)
    cat(paste0(rep("_",50), collapse = ""), file = "./results/used_parameters.txt",append = T, sep = "\n")
  }
  if (impute == "yes"){
    print("Missing values will be imputed using mice. A maximum proportion of ")
    print(paste0("A maximum proportion of ",NAs_prop," is accepted"))
    cols <- c("protein_group","exp_group")
    dataset <- dataset %>%
      group_by(across(all_of(cols))) %>%
      mutate(count_NA = sum(is.na(normalized_intensity))) %>%
      mutate(per_NA = sum(is.na(normalized_intensity)) /sum(sum(is.na(normalized_intensity)) + sum(!is.na(normalized_intensity)))) %>%
      ungroup() %>%
      mutate(was_NA = ifelse(is.na(normalized_intensity) == T,"Was NA","Was not NA"))
    
    # Filter pg with less than NAs_prop
    dataset <- dataset %>%
      filter(per_NA <= NAs_prop)
    
    # Reshape the data to impute
    to_imput <- dcast(dataset, 
                      protein_group ~ sample_name, value.var="normalized_intensity")
    
    rownames(to_imput) <- to_imput$protein_group
    to_imput <- to_imput[,-1]
    
    # Impute the data
    imputed_data <- complete(mice(to_imput, m = 5, method = "pmm", seed = 500))
    
    # Reshape the data again to long format
    samples <- rep(colnames(imputed_data), each = length(unique(rownames(imputed_data))))
    pgs <- rep(rownames(imputed_data), times = length(unique(samples)))
    values <- unlist(as.vector(imputed_data), use.names = F)
    
    imputed_data_long <- data_frame("sample_name" = samples,
                                    "protein_group" = pgs,
                                    "provisional" = values)
    
    dataset[,"normalized_intensity_before_imputation"] <- dataset[,"normalized_intensity"]
    dataset <- subset(dataset, select = -c(normalized_intensity))
    
    imputed_data_long[,"normalized_intensity"] <- imputed_data_long[,"provisional"]
    imputed_data_long <- subset(imputed_data_long, select = -c(provisional))
    
    imputed_data <- merge(imputed_data_long, dataset, by = c("sample_name","protein_group"))
    imputed_data <- imputed_data %>%
      relocate(normalized_intensity_before_imputation, .after = intens)
    
    cols <- c("protein_group")
    imputed_data <- imputed_data %>%
      group_by(across(all_of(cols))) %>%
      filter(n() >= n) %>%
      ungroup()
    cat("tim - for NA imputation",file = "./results/used_parameters.txt",sep = "\n", append = T)
    cat(paste0("impute == yes ",unique(NAs_prop)), file = "./results/used_parameters.txt", sep = "\n", append = T)
    cat(paste0(rep("_",50), collapse = ""), file = "./results/used_parameters.txt",append = T, sep = "\n")
  }
  return(imputed_data)
}

pair_wise_pca <- function(dataset,group1,group2,column,to_colour){
  print(paste0("Plotting ",group1," and ",group2," in the PCA"))
  part1 <- dataset[dataset$exp_group == group1,]
  part2 <- dataset[dataset$exp_group == group2,]
  data <- rbind(part1,part2)
  data <- dcast(data, protein_group ~ 
                  sample_name,value.var="normalized_intensity", fun.aggregate = median)
  data <- t(data[,c(2:ncol(data))])
  pca <- prcomp(data, scale. = TRUE, center = TRUE)
  samples_afer_batch <- intersect(rownames(data), to_colour$sample_name)
  to_colour2 <- subset(to_colour, sample_name %in% c(samples_afer_batch))
  pca2_graph <- autoplot(pca, data = to_colour2, colour = "exp_group",
                         frame = T)+
    scale_fill_manual(values = cbp1) +
    scale_color_manual(values = rep("black",9))
  return(pca2_graph)}

make_all_contrasts <- function(design){
  group <- unique(as.character(colnames(design)))
  cb   <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  contrasts<- limma::makeContrasts(contrasts=cb, levels=group)
  colnames(contrasts) <- gsub("-", "_vs_", colnames(contrasts))
  return(contrasts)}

upsidedown <- function(contmat, comparisons_to_change=NULL){
  if (is.null(comparisons_to_change)){
    contmat <- contmat}
  if (!is.null(comparisons_to_change)){
    for (i in 1:length(contmat)){
      if (colnames(contmat)[i] %in% comparisons_to_change){
        contmat[,i] <- (contmat[,i] * (-1))
      }
    }
    for (j in 1:length(colnames(contmat))){
      if (colnames(contmat)[j] %in%  comparisons_to_change){
        before <- str_split(colnames(contmat)[j], "_vs_", simplify = T)[1]
        after <- str_split(colnames(contmat)[j], "_vs_", simplify = T)[2]
        colnames(contmat)[j] <- paste0(after,"_vs_",before,collapse = "")
      }
    }
  }
  return(contmat)
}

tt_extractor <- function(fit1,annotation){
  coefficients <- colnames(head(fit1@.Data)[[1]])
  tt_list <- list()
  for (i in 1:length(colnames(head(fit1@.Data)[[1]]))){
    top_table_unique <- topTable(fit = fit1, coef = coefficients[[i]], number = Inf)
    top_table_unique <- merge(top_table_unique, annotation, by = "row.names")
    colnames(top_table_unique)[1] <- "protein_group"
    top_table_unique <- top_table_unique[, !duplicated(colnames(top_table_unique))]
    top_table_unique <- relocate(.data = top_table_unique, c(accession), .after = protein_group)
    tt_list[[i]] <- top_table_unique
  }
  names(tt_list) <- coefficients
  names(tt_list) <- gsub(pattern = " ", replacement = "_", x = names(tt_list))
  return(tt_list)
}

tt_list_cleaner <- function(list, meta_data){
  final_names <- names(list)
  final_list <- list()
  for (i in 1:length(list)){
    cleanedtt <- as.data.frame(list[[i]])
    cleanedtt <- subset(cleanedtt, select = -c(AveExpr,t,B))
    cleanedtt$Fold_Change <- 2^cleanedtt$logFC
    cleanedtt <- relocate(.data = cleanedtt, Fold_Change, .after = logFC)
    
    # Get the experimental group names
    both_samp <- final_names[i]
    sample_1 <- str_split(both_samp, "_vs_", simplify = T)[1]
    sample_2 <- str_split(both_samp, "_vs_", simplify = T)[2]
    
    # Corresponding to meta_data 2
    meta_data2 <- meta_data[meta_data$sample_name %in% to_colour$sample_name,]
    
    # Correspondance
    
    correspond1 <- meta_data2$sample_name[meta_data2$exp_group == sample_1]
    correspond2 <- meta_data2$sample_name[meta_data2$exp_group == sample_2]
    correspond <- c(correspond1,correspond2)
    
    # Corresponding to meta_data 2
    meta_data2 <- meta_data2[meta_data2$sample_name %in% correspond,]
    
    # Get the data
    to_bind <- as.data.frame(t(maxquant_mat))
    rownames(to_bind) <- maxquant_imp_to_PCA$protein_group
    to_bind <- to_bind[colnames(to_bind) %in% c(correspond1,correspond2)]
    to_bind <- to_bind %>% 
      select(all_of(correspond))
    to_bind$protein_group <- rownames(to_bind)
    
    # final merge
    cleanedtt <- merge(to_bind, cleanedtt, by = "protein_group")
    final_list[[i]] <- cleanedtt
  }
  names(final_list) <- final_names
  return(final_list)
}
#----

### Work with the dataset ----
# Retrieve only intensity columns form maxquant dataset
maxquant <- use_this_intensity(dataset = maxquant, pattern = "LFQ intensity")

# Save sample names
sample_names <- colnames(maxquant[-c(1,2)])

# Transform 0s to NAs
maxquant <- zero_to_NA(dataset = maxquant, patterns = sample_names)

# Remove contaminants
maxquant_clean <- removecontaminants(dataset = maxquant, contaminants = cont, accession_name = "accession")

# transform to log2
maxquant_clean <- log2_to_pattern(dataset = maxquant_clean, patterns = sample_names)
#----

### Data Quality----
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
          "#00B159", "#FCD612", "#FF2F03", "#03D3FF",
          "#2FF923", "#FF03B4", "#f08a0c", "#0d407f",
          "#037A68", "#510840", "#F70D1A", "#e0218a")


# Data to long format
long_format <- maxquant_clean %>%
  pivot_longer(cols = all_of(sample_names),
               names_to = "sample_name",
               values_to = "intens")
# Add metadat to long format
long_format <- merge(long_format, meta_data, by = "sample_name")

# Add median calculation (by sample)
long_format_median <- long_format %>%
  group_by(sample_name) %>%
  summarise(MED = median(intens, na.rm=T))

long_format <- merge(x = long_format, y = long_format_median,
                     by = "sample_name")

# Boxplot
intensity_boxplots <- ggplot(long_format, mapping = aes(x = sample_name, y = intens))+
  geom_boxplot(data = long_format, mapping = aes(x = sample_name, y = intens, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots

# Completenes
completeness <- long_format %>%
  group_by(sample_name) %>%
  summarize(count_na = 100 - (((sum(is.na(intens))/nrow(maxquant_clean))))*100)
completeness <- merge(completeness, meta_data, by = "sample_name")
completeness <- completeness %>%
  mutate(plex = 
           substr(sample_name, 1,5))

completeness_barplot <- ggplot(completeness, mapping = aes(y = count_na, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = plex))+
  geom_hline(data = completeness ,mapping = aes(yintercept = mean(count_na), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("% of completeness")+
  ggtitle("Completeness of the data per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
completeness_barplot

# Number of proteins per sample
nprot <- long_format %>%
  group_by(sample_name) %>%
  summarize(count_prot = sum(!is.na(intens)))
nprot <- merge(nprot, meta_data, by = "sample_name")
nprot <- nprot %>%
  mutate(plex = 
           substr(sample_name, 1,5))

number_of_prot <- ggplot(nprot, mapping = aes(y = count_prot, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = plex))+
  geom_hline(data = nprot ,mapping = aes(yintercept = mean(count_prot), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("# of proteins")+
  ggtitle("Proteins per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
number_of_prot

# NA per protein
naprot <- long_format %>%
  group_by(protein_group) %>%
  summarize(count_prot = sum(is.na(intens)))

# NAs density per protein
na_density <- ggplot(data = naprot, mapping = aes(x = count_prot))+
  geom_histogram(binwidth = 1, bins = 1)+
  theme_bw()+
  ggtitle("NAs density per protein groups")+
  xlab("# NAs")+
  ylab("# proteins")+
  theme(legend.position= "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
na_density
#----

### Remove samples----
long_format <- remove_samp(dataset = long_format, samples = c("cec23","cec7","ns14"))
#----

### Normalization by median----
maxquant_median <- long_format

# Get global median
median_all <- maxquant_median$intens

# Do the calculation
maxquant_median$normalized_intensity <- (maxquant_median$intens - maxquant_median$MED) + median_all 

intensity_boxplots_norm <- ggplot(maxquant_median, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = maxquant_median, mapping = aes(x = sample_name, y = normalized_intensity, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

intensity_boxplots_norm
#----

### NA imputation----
maxquant_clean_median_imput <- tim(impute = "yes", dataset = maxquant_median, NAs_prop = 0.5)
#----

### PCA----
# Change long to wide format
maxquant_imp_to_PCA <- dcast(maxquant_clean_median_imput, 
                             protein_group
                             ~ sample_name,value.var="normalized_intensity",
                             fun.aggregate = median)

maxquant_mat <- t(maxquant_imp_to_PCA[,c(2:ncol(maxquant_imp_to_PCA))])

# PCA
pca1 <- prcomp(maxquant_mat, scale. = TRUE, center = TRUE)

# PCA colouring
to_colour <- as.data.frame(maxquant_mat)
to_colour$sample_name <- rownames(to_colour)
  
# Add metadata to the to_colour dataset 
to_colour <- merge(to_colour, meta_data, by = "sample_name")

# Plot the graph
pca1_graph <- autoplot(pca1, data = to_colour, colour = "exp_group",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca1_graph
#----

### Pair-wise PCA----
# Check for the experimental groups
unique(to_colour$exp_group)

# Plot PCAs
pair_wise_pca(dataset = maxquant_clean_median_imput, group1 = "Squamous",
              group2 = "Basal", to_colour = to_colour)
#----

### limma---- 
# Retrive expression matrix
expression_matrix <- as.data.frame((dcast(maxquant_clean_median_imput, 
                                          protein_group ~ sample_name,value.var="normalized_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)

# Create design and contrast matrix
to_colour <- subset(to_colour, sample_name %in% colnames(expression_matrix))
to_colour <- to_colour[match(colnames(expression_matrix), to_colour$sample_name)]
to_colour <- subset(to_colour, select = sample_name)
to_colour <- merge(to_colour, meta_data)

groups <- to_colour$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

# Annotation
annotation <- subset(maxquant, select = c(protein_group,accession))
annotation <- distinct(annotation)
rownames(annotation) <- annotation$protein_group

# Fit the model
fit <- lmFit(expression_matrix, design = design)
#----

### limma contrasts----
# Prepare all possile contrasts
contrasts_all <- make_all_contrasts(design)

# check the direction
colnames(contrasts_all)

# Reverse the desired contrasts
contrasts_all <- upsidedown(contmat = contrasts_all, comparisons_to_change = "Normal_Sample_vs_Squamous")
contrasts_all
#----

### Fit contrasts to the model----
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)
#----

### Extract TT----
# For all contrasts
tt <- topTable(fit1, number = Inf)
tt <- merge(tt, annotation, by = "row.names")
colnames(tt)[1] <- "protein_group"
tt <- tt[, !duplicated(colnames(tt))]
tt <- tt %>%
  relocate(accession, .after = protein_group)
write.table(tt, "./results/TopTable_all.tsv", row.names = F, sep = "\t", dec = ".")

# Extract all Top Tables
tt_all <- tt_extractor(fit1 = fit1, annotation = annotation)
#----

### Clear TTs----
tt_all <- tt_list_cleaner(tt_all, meta_data = meta_data)
#----

### Write xlsx file----
write.xlsx(tt_all, file = "./results/final_toptables.xlsx")
#----

### Heatmap----
# Filter the data
ttplot <- tt
ttplot <- ttplot %>% 
  filter(adj.P.Val < 0.05)

# Create rownnames matrix for the plot
heat_matrix <- as.matrix(ttplot[,grep("vs",colnames(ttplot))])
rownames(heat_matrix) <- ttplot$accession

# Plot the heatmap
all_heatmap <- heatmap.2(x = heat_matrix,
                         trace = "none", density.info = "none",
                         main = "Differential protein expression", scale = "row", )

#----



