# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings", version = "3.15", force = TRUE)

install.packages(c('bioseq', 'qpcR', 'rgl', 'glmnet', 'e1071', 'tidyverse', 'ggplot2', 'reshape2', 'data.table', 'randomForest'))

library(Biostrings)
library(rentrez)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(bioseq)
library(randomForest)
library(qpcR)
library(rgl)
library(e1071)

# Set the working directory to the project root directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..") # Move to the project root directory

# Searches the terms that can be used to extract genes with regards to the database terms
entrez_db_searchable("nuccore")

# Extracting HLA-A and B from the database
HLA_A_gene <- entrez_search(db = "nuccore", term = "Human[ORGN] AND HLA-A[GENE]", retmax = 1000, use_history = TRUE)
HLA_B_gene <- entrez_search(db = "nuccore", term = "Human[ORGN] AND HLA-B[GENE]", retmax = 1000, use_history = TRUE)

# DATA Acquisition for HLA-A sequences
web_history_HLA_A_seq <- entrez_fetch(db = "nuccore", web_history = HLA_A_gene$web_history, rettype = "fasta")
write(web_history_HLA_A_seq, "data/HLA_A_seq.fasta", sep = "\n")
HLA_A_seq_stingset <- readDNAStringSet("data/HLA_A_seq.fasta")

# Creating a dataframe for the HLA-A
df_HLA_A_seq <- data.frame(HLA_A_Title = names(HLA_A_seq_stingset), HLA_A_Sequence = paste(HLA_A_seq_stingset))

# Filter and clean the HLA-A sequences
df_HLA_A_seq <- df_HLA_A_seq[nchar(df_HLA_A_seq$HLA_A_Sequence) > 500, ]
df_HLA_A_seq <- df_HLA_A_seq[grepl("[^N]*", df_HLA_A_seq$HLA_A_Sequence), , drop = FALSE]

# DATA Acquisition for HLA-B sequences
web_history_HLA_B_seq <- entrez_fetch(db = "nuccore", web_history = HLA_B_gene$web_history, rettype = "fasta")
write(web_history_HLA_B_seq, "data/HLA_B_seq.fasta", sep = "\n")
HLA_B_seq_stingset <- readDNAStringSet("data/HLA_B_seq.fasta")

# Creating a dataframe for the HLA-B
df_HLA_B_seq <- data.frame(HLA_B_Title = names(HLA_B_seq_stingset), HLA_B_Sequence = paste(HLA_B_seq_stingset))

# Filter and clean the HLA-B sequences
df_HLA_B_seq <- df_HLA_B_seq[nchar(df_HLA_B_seq$HLA_B_Sequence) > 500, ]
df_HLA_B_seq <- df_HLA_B_seq[grepl("[^N]*", df_HLA_B_seq$HLA_B_Sequence), , drop = FALSE]

# Calculate nucleotide frequencies and generate visualizations
df_GC_content <- qpcR:::cbind.na(
  data.frame(x = letterFrequency(DNAStringSet(df_HLA_B_seq$HLA_B_Sequence), as.prob = TRUE, letters = c("GC"))),
  data.frame(x = letterFrequency(DNAStringSet(df_HLA_A_seq$HLA_A_Sequence), as.prob = TRUE, letters = c("GC")))
)

df_AT_content <- qpcR:::cbind.na(
  data.frame(x = letterFrequency(DNAStringSet(df_HLA_B_seq$HLA_B_Sequence), as.prob = TRUE, letters = c("AT"))),
  data.frame(x = letterFrequency(DNAStringSet(df_HLA_A_seq$HLA_A_Sequence), as.prob = TRUE, letters = c("AT")))
)

# Creating the Boxplot for nucleotide frequencies
df_GC_content <- cbind(Dinucleotide = 'GC', df_GC_content)
df_AT_content <- cbind(Dinucleotide = 'AT', df_AT_content)
DF_GC_AT_content1 <- rbind(df_GC_content, df_AT_content)

df_long <- reshape2::melt(DF_GC_AT_content1, id.vars = c("Dinucleotide"), variable.name = "Gene", value.name = "Percentage")
df1_long <- na.omit(df_long)

Frequency_Boxplot <- ggplot(df1_long, aes(x = Gene, y = Percentage, fill = Gene)) +
  geom_boxplot() + labs(title = "AT & GC Frequency") + facet_wrap(~Dinucleotide) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the boxplot to the output directory
ggsave("output/Frequency_Boxplot.png", plot = Frequency_Boxplot)

# Mean length barplot for HLA-A and B
HLA_A_Mean_length <- mean(width(HLA_A_seq_stingset))
HLA_B_Mean_length <- mean(width(HLA_B_seq_stingset))
Mean_len <- data.frame(Gene = c("HLA-A", "HLA-B"), Sequence_Mean_length = c(HLA_A_Mean_length, HLA_B_Mean_length))

Mean_length_barplot <- ggplot(data = Mean_len, aes(x = Gene, y = Sequence_Mean_length)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  labs(y = "Sequence Mean Length (bp)", x = "Gene", title = "Average Length of HLA-A and HLA-B Sequences") +
  geom_text(aes(label = Sequence_Mean_length), vjust = 1.6, color = "white", size = 3.5) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the mean length barplot to the output directory
ggsave("output/Mean_length_barplot.png", plot = Mean_length_barplot)

# Prepare data for machine learning classification
df_HLA_B_seq$Gene <- 'HLA-B'
df_HLA_A_seq$Gene <- 'HLA-A'
df_Combined <- data.table(merge(df_HLA_B_seq, df_HLA_A_seq, by = "Gene", all = TRUE))
df_final <- melt(df_Combined, id.vars = "Gene", measure.vars = c("HLA_A_Sequence", "HLA_B_Sequence"), value.name = "Sequence")
df_final <- na.omit(df_final)
df_final$Sequence <- DNAStringSet(df_final$Sequence)
df_final <- cbind(df_final, as.data.frame(letterFrequency(df_final$Sequence, letters = c("A", "C", "G", "T"))))
df_final <- cbind(df_final, as.data.frame(dinucleotideFrequency(df_final$Sequence, as.prob = TRUE)))
df_final <- cbind(df_final, as.data.frame(trinucleotideFrequency(df_final$Sequence, as.prob = TRUE)))

# Split data into training and validation sets
set.seed(200)
dfValidation <- df_final %>%
  group_by(Gene) %>%
  sample_n(1000)

set.seed(13)
df_Training <- df_final %>%
  filter(!Sequence %in% dfValidation$Sequence) %>%
  group_by(Gene) %>%
  sample_n(61)

# Train Random Forest Classifier
gene_classifier <- randomForest::randomForest(x = df_Training[, 3:89], y = as.factor(df_Training$Gene), ntree = 50, importance = TRUE)

# Evaluate using Random Forest
predictValidation <- predict(gene_classifier, dfValidation[, 3:89])
confusion_matrix <- table(observed = dfValidation$Gene, predicted = predictValidation)
CM_RF_table <- data.frame(confusion_matrix)

# Plot Random Forest confusion matrix
RF_plot <- ggplot(data = CM_RF_table, mapping = aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "Random Forest") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the Random Forest plot
ggsave("output/Random_Forest_Plot.png", plot = RF_plot)

# Train Naive Bayes Classifier
gene_classifier_Bayes <- naiveBayes(df_Training[, 3:89], df_Training$Gene)

# Evaluate using Naive Bayes
predictValidation <- predict(gene_classifier_Bayes, newdata = dfValidation[, 3:89])
cm_Bayes <- table(dfValidation$Gene, predictValidation)
cm_Bayes_table <- data.frame(cm_Bayes)

# Plot Naive Bayes confusion matrix
Bayes_plot <- ggplot(data = cm_Bayes_table, mapping = aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "Naive Bayes") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the Naive Bayes plot
ggsave("output/Naive_Bayes_Plot.png", plot = Bayes_plot)
