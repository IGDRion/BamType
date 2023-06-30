# Set the directory path where the CSV files are located
directory <- "/groups/dog/stage/victor/bamslam/csv/subsample"

# Get the list of file names with the ".csv" extension
file_names <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty dataframe
merged_df <- data.frame()

# Iterate over the file names and read each CSV file
for (file in file_names) {
  if (file != paste0(directory,"/all_subsample.csv")){
    df <- read.csv(file)
    merged_df <- rbind(merged_df, df)
  }
}

# Add column of condition (sensitive / resistant)
merged_df$condition <- ifelse(grepl("S_R1|S1_R1|S2_R1|S3_R1", merged_df$sample), "sensitive",
                         ifelse(grepl("R_R1|R1_R1|R2_R1|R3_R1", merged_df$sample), "resistant", NA))

# Add column of cancer
merged_df$cancer <- ifelse(grepl("^501Mel", merged_df$sample), "Melanoma",
                        ifelse(grepl("^PC3", merged_df$sample), "Prostate",
                               ifelse(grepl("^U251", merged_df$sample), "Glioblastoma",
                                      ifelse(grepl("^ADCA72", merged_df$sample), "Lung", NA))))

# Keep only rows that are PCG and lncRNA
merged_df <- subset(merged_df, transcript_biotype == "protein_coding" | transcript_biotype == "lncRNA")

merged_df$sample <- as.factor(merged_df$sample)
merged_df$condition <- as.factor(merged_df$condition)
merged_df$cancer <- as.factor(merged_df$cancer)
merged_df$transcript_biotype <- as.factor(merged_df$transcript_biotype)

# Print the merged dataframe
head(merged_df)

# Write merged csv file
write.table(merged_df, file = paste0(directory,"/all_sample.csv"), sep=",", quote=F, col.names = TRUE, row.names = FALSE)
