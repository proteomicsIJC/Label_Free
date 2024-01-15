##########################
### Label Free Worflow ###
##########################

### Load libraries, WD and functions
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
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

# Functions deffinition
## General proteomics functions
source("../functions/general/zero_to_na_label_free.R")
source("../functions/general/presence_vs_no_presence.R")
source("../functions/TMT_MaxQuant/proteinGroupsCleaner.R")
source("../functions/general/log2_to_pattern_label_free.R")
source("../functions/general/remove_samp.R")
source("../functions/general/tim.R")
source("../functions/general/make_all_contrasts.R")
source("../functions/general/upsidedown.R")

## Label Free MaxQuant proteomics
source("../functions/Label_Free_MaxQuant/use_this_intensity.R")

## Folder system
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))

### Data importation
raw_data <- readr::read_tsv("./raw_data/proteinGroups.txt")

# Meta data
meta_data <- openxlsx::read.xlsx("./raw_data/meta_data.xlsx", sheet = 1)

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

### Data work
# Remove Contaminants
raw_data <- raw_data %>% 
  filter(is.na(`Potential contaminant`),
         is.na(Reverse),
         is.na(`Only identified by site`))
maxquant_clean <- raw_data

# Retrieve only intensity columns form maxquant dataset
maxquant_clean <- use_this_intensity(dataset = maxquant_clean, pattern = "LFQ intensity")
for (i in 1:length(colnames(maxquant_clean))){
  colnames(maxquant_clean)[i] <- gsub(pattern = "LFQ intensity ", replacement = "", colnames(maxquant_clean)[i])
}
colnames(maxquant_clean)[1] <- "protein_group"

# Save sample names
sample_names <- colnames(maxquant_clean[-c(1)])

# Transform 0s to NAs
maxquant <- zero_to_NA_label_free(dataset = maxquant_clean, patterns = sample_names)

# transform to log2
maxquant_clean <- log2_to_pattern_label_free(dataset = maxquant_clean, patterns = sample_names)


### Data Quality
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

### Remove samples
long_format <- remove_samp(dataset = long_format, samples = c("CEC23","NS7","NS14"))

### Normalization by median
maxquant_median <- long_format

# Get global median
median_all <- median(maxquant_median$intens, na.rm = T)

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

### NA imputation
## create exp_group
maxquant_median$exp_group <- "group1"

## Do the imputation
maxquant_clean_median_imput <- tim(impute = "yes", dataset = maxquant_median, NAs_prop = 0.5)

#### HERE DO THE IMUTATION WITH THE MATRIX AND SEE IF IT WORKS FINE !!!!!!!
imp <- impute.knn(as.matrix(mydata.ff.50[, 142:165]))$data
mydata.ff.50[, 142:165] <- imp
#####

### PCA
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
pca1_graph <- autoplot(pca1, data = to_colour, colour = "sample_group",
                       frame = T)+
  scale_fill_manual(values = cbp1) +
  scale_color_manual(values = rep("black",9))
pca1_graph


### limma
# Retrive expression matrix
expression_matrix <- as.data.frame((dcast(maxquant_clean_median_imput, 
                                          protein_group ~ sample_name,value.var="normalized_intensity", fun.aggregate = median)))

rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[,-1]
expression_matrix <- expression_matrix %>% 
  mutate_if(is.character, as.numeric)

## Create design and contrast matrix HERE PUT AN OBJECT CALLED META DATA TRACKER !!!!!!!!! 
to_colour <- subset(to_colour, sample_name %in% colnames(expression_matrix))
to_colour <- to_colour[match(colnames(expression_matrix), to_colour$sample_name)]
to_colour <- subset(to_colour, select = sample_name)
to_colour <- merge(to_colour, meta_data)

groups <- to_colour$exp_group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

## Annotation
annotation <- subset(maxquant, select = c(protein_group,accession))
annotation <- distinct(annotation)
rownames(annotation) <- annotation$protein_group

## Fit the model
fit <- lmFit(expression_matrix, design = design)

# limma contrasts
## Prepare all possile contrasts
contrasts_all <- make_all_contrasts(design)

## check the direction
colnames(contrasts_all)

## Reverse the desired contrasts
contrasts_all <- upsidedown(contmat = contrasts_all, comparisons_to_change = "Normal_Sample_vs_Squamous")
contrasts_all

## Fit contrasts to the model
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)

# Extract TTs and write xlsx
## For all contrasts
tt <- topTable(fit1, number = Inf)
tt <- merge(tt, annotation, by = "row.names")
colnames(tt)[1] <- "protein_group"
tt <- tt[, !duplicated(colnames(tt))]

## Write xlsx

tt <- tt %>%
  relocate(accession, .after = protein_group)
write.table(tt, "./results/TopTable_all.tsv", row.names = F, sep = "\t", dec = ".")

