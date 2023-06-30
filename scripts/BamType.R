main <- function() {
  
    # ---------------------------------------------------------------------------- #
    #                      Loading packages & data importation                     #
    # ---------------------------------------------------------------------------- #

    args = commandArgs(trailingOnly=TRUE)
    type <- strsplit(grep('--type*', args, value = TRUE), split = '=')[[1]][[2]] # cdna or rna
    bamfile <- strsplit(grep('--bam*', args, value = TRUE), split = '=')[[1]][[2]] # alignment file
    output <- strsplit(grep('--output*', args, value = TRUE), split = '=')[[1]][[2]] # path of output
    gtffile <- strsplit(grep('--gtf*', args, value = TRUE), split = '=')[[1]][[2]] # extended annotation
    seqsummaryfile <- read.table(file=(strsplit(grep('--seqsum*', args, value = TRUE), split = '=')[[1]][[2]]), header=TRUE) # sequencing summary file (contains qscore information)

    sample_name <- basename(bamfile)
    sample_name <- gsub("\\.bam$", "", sample_name)

    suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(dplyr, warn.conflicts = FALSE)
    library(ggplot2)
    library(hexbin)
    library(gridExtra)
    library(rtracklayer)
    })

    options(dplyr.summarise.inform = FALSE)

    message("Starting sample: ", bamfile)

    
    # Import bam
    bam <- readGAlignments(bamfile, use.names = TRUE,
                            param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                                what = c("qname","flag","mapq")))
    message("Imported bam file")




    # ---------------------------------------------------------------------------- #
    #                             Modifying BAM object                             #
    # ---------------------------------------------------------------------------- #


    # Expand CIGAR strings
    cigaropts = cigarOpTable(bam@cigar)

    for (col in colnames(cigaropts)) {
    mcols(bam)[[paste0("nbr", col)]] = cigaropts[, col]
    }

    message("Expanded CIGAR strings")



    # Transcript origin information
    bam_data <- as.data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>% # change the GAlignment object to a data frame object
    dplyr::select(-cigar) %>%
    dplyr::select(-njunc) %>%
    dplyr::mutate(novel_known = if_else(grepl('tx.', seqnames) == T, 'novel', 'known')) # create a column 'novel_known', and for each transcript if it novel (start with 'tx.'), 'novel_known' takes value 'novel', else it takes 'known'

    bam_data$novel_known <- as.factor(bam_data$novel_known)

    nbr_align <- length(bam_data$qname)
    nbr_read <- length(unique(bam_data$qname))

    message("Transcript status column added")


    # Creation of a length column (length of reference transcript)
    lengths <- as.data.frame(bam@seqinfo) %>% 
    dplyr::select(-isCircular) %>% 
    dplyr::select(-genome) %>% 
    tibble::rownames_to_column("reference")

    bam_data <- merge(bam_data, lengths, by.x="seqnames", by.y="reference", all.x=TRUE) # Add seqlength column on bam_data

    rm(bam) # Clear memory

    message("Transcript length column added")


    # Merging the bam and the sequencing summary file to add qscore information for each read (merge by read id (qname))
    bam_data <- merge(x = bam_data, y = seqsummaryfile[ ,c("read_id","mean_qscore_template")], all = TRUE, by.x = "qname", by.y = "read_id")

    message("Qscore column created")

    rm(seqsummaryfile) # Clear memory



    # Transcript biotype information
    # import is from rtracklayer package. Keeping only type=="transcript" rows ; and transcript_id + transcript_biotype columns for merging
    gtf <- import(gtffile)
    gtf <- as.data.frame(gtf) %>%
    subset(type=="transcript") %>%
    dplyr::select(transcript_id, transcript_biotype)

    # Merging the bam and the GTF to add transcript biotype information (merge by transcript id)
    bam_data <- merge(x = bam_data, y = gtf, all = TRUE, by.x = "seqnames", by.y = "transcript_id")
    bam_data$transcript_biotype <- as.factor(bam_data$transcript_biotype)

    rm(gtf) # Clear memory
    nbr_transcript_unmapped <- sum(as.data.frame(is.na(bam_data$novel_known)))

    message("Biotype column created")



    # Select alignements with their flagstats (also gets rid of NA rows)
    if (type == "cdna") {
    bam_data <- subset(bam_data, flag == 0 | flag == 16 | flag == 256 | flag == 272)
    } else if (type == "rna") {
    bam_data <- subset(bam_data, flag == 0 | flag == 256)
    } else {
    print("Sequencing type missing. Please enter either: cdna rna")
    }

    new_nbr_align <- length(bam_data$qname)
    new_nbr_read <- length(unique(bam_data$qname))

    message("Flagstat filtering done")



    # Creation of new informative columns with CIGAR string
    bam_data <- bam_data %>% 
    dplyr::mutate(alignedLength = nbrM + nbrI) %>% 
    dplyr::mutate(readLength = nbrS + nbrH + nbrM + nbrI) %>% 
    dplyr::mutate(alignedFraction = alignedLength/readLength) %>% #This is the coverage of the read (how much of the read is aligned to the transcript)
    dplyr::mutate(accuracy=(nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD))

    
    # Creation of transcript coverage column (how much of the transcript is covered)
    bam_data <- bam_data %>% 
    dplyr::mutate(coverage=width/seqlengths)

    bam_data$coverage <- as.numeric(bam_data$coverage)

    message("Calculated transcripts and reads coverages")


    # Creation of number of secondary alignments column
    alignments <- bam_data %>% 
    dplyr::group_by(qname) %>% 
    dplyr::summarise(nbrSecondary = n()-1) %>% 
    dplyr::rename(read = qname)

    bam_data <- merge(bam_data, alignments, by.x="qname", by.y="read", all.x=TRUE)

    alignments <- alignments %>% 
    dplyr::group_by(nbrSecondary) %>% 
    dplyr::summarise(total = n()) %>%
    dplyr::mutate(prop = total / sum(total))

    bam_data$nbrSecondary <- as.factor(bam_data$nbrSecondary)
    alignments$nbrSecondary <- as.factor(alignments$nbrSecondary)

    message("Calculated number of secondary alignments")

    # Export whole data file

    if (!dir.exists(dirname(output))){
                stop(paste0(dirname(output))," directory does not exist! Exiting the script")}
        else{
            if (!dir.exists(output)){
                dir.create(output)}
        }

    if (!dir.exists(paste0(output,"csv"))){
                dir.create(paste0(output,"csv"))}

    if (!dir.exists(paste0(output,"csv/all_data"))){
                dir.create(paste0(output,"csv/all_data"))}
    
    bam_export <- subset(bam_data, select=c("qname", "seqnames", "start", "end", "flag", "mapq", "AS", "novel_known", "transcript_biotype", "alignedLength", "readLength", "alignedFraction", "accuracy", "seqlengths", "coverage", "nbrSecondary"))
    write.csv(bam_export, file = paste0(output,"csv/all_data/",sample_name,"_data.csv"), sep=",", quote=F, col.names = T, row.names=F)

    message("Exported whole data file")




    # ---------------------------------------------------------------------------- #
    #    Making plots for subgroups + exporting primary alignments data & stats    #
    # ---------------------------------------------------------------------------- #


    #Plot subgroups
    make_plots <- function(){
    
    subgroups_order <- c()

    cov_med_vec <- c()
    read_len_med_vec <- c()
    prop_full_tr_vec <- c()
    align_frac_med_vec <- c()
    accuracy_med_vec <- c()
    sec_align_mean_vec <- c()
    qscore_med_vec <- c()

    plots1 <- list()
    plots2 <- list()
    plots3 <- list()
    plots4 <- list()
    plots5 <- list()

    status_vector = c("known", "novel", "full")
    biotype_vector = c("protein_coding", "lncRNA")
    count <- 1

    # Create filtered bam data frame
    for (i in 1:length(status_vector)){
        status <- status_vector[i]

        for (j in 1:length(biotype_vector)){
            biotype <- biotype_vector[j]

            if (status == "full"){
                bam_data_filtered <- bam_data
            }
            else {
                bam_data_filtered <- bam_data %>%
                dplyr::filter(grepl(status, novel_known)) %>%
                dplyr::filter(grepl(biotype, transcript_biotype))
            }

            if (status=="full") {
                message("Starting plot generation for full data")
            }
            else {
                message(paste0("Starting plot generation for subset ",status,"_",biotype))
            }
            
            # Select the best AS per read, but prioritise original primary alignment if present
            bam_primary <- bam_data_filtered %>% 
            dplyr::group_by(qname) %>% 
            dplyr::arrange(tp) %>% 
            dplyr::arrange(qname, desc(AS)) %>% 
            dplyr::slice(1)

            bam_primary <- bam_primary %>% 
            dplyr::mutate(above=coverage>0.95)


            # Make alignement data frame for the subgroup
            alignments <- dplyr::select(bam_primary, nbrSecondary, qname)

            alignments <- alignments %>%
            dplyr::group_by(nbrSecondary) %>% 
            dplyr::summarise(total = n()) %>%
            dplyr::mutate(prop = total / sum(total))

            alignments$nbrSecondary <- as.factor(alignments$nbrSecondary)


            # For each transcript associate the highest coverage
            bam_per_unique_transcript <- bam_primary %>% 
            dplyr::group_by(seqnames) %>% 
            summarise(coverage = median(coverage, na.rm = TRUE))


            # For each transcript associate its length
            length_per_unique_transcript <- bam_primary%>% 
            dplyr::group_by(seqnames) %>% 
            summarise(seqlengths = max(seqlengths))

            message("bam_primary created")


            # Subgroups stats for common plots

            if(status=="full"){
                # Proportion information and number of read
                nbr_read_known_pcg <- nrow(as.data.frame(bam_primary[bam_primary$transcript_biotype=="protein_coding" & bam_primary$novel_known=="known", ]))
                nbr_read_novel_pcg <- nrow(as.data.frame(bam_primary[bam_primary$transcript_biotype=="protein_coding" & bam_primary$novel_known=="novel", ]))
                nbr_read_known_lnc <- nrow(as.data.frame(bam_primary[bam_primary$transcript_biotype=="lncRNA" & bam_primary$novel_known=="known", ]))
                nbr_read_novel_lnc <- nrow(as.data.frame(bam_primary[bam_primary$transcript_biotype=="lncRNA" & bam_primary$novel_known=="novel", ]))
                nbr_read_known_oth <- nrow(as.data.frame(bam_primary[(bam_primary$transcript_biotype!="protein_coding" & bam_primary$transcript_biotype!="lncRNA") & bam_primary$novel_known=="known", ]))
                nbr_read_novel_oth <- nrow(as.data.frame(bam_primary[(bam_primary$transcript_biotype!="protein_coding" & bam_primary$transcript_biotype!="lncRNA") & bam_primary$novel_known=="novel", ]))

                prop_read_known_pcg <- nbr_read_known_pcg/new_nbr_read
                prop_read_novel_pcg <- nbr_read_novel_pcg/new_nbr_read
                prop_read_known_lnc <- nbr_read_known_lnc/new_nbr_read
                prop_read_novel_lnc <- nbr_read_novel_lnc/new_nbr_read
                prop_read_known_oth <- nbr_read_known_oth/new_nbr_read
                prop_read_novel_oth <- nbr_read_novel_oth/new_nbr_read

                nbr_unique_tr<- nrow(bam_per_unique_transcript)
                med_cov_unique_tr <- median(bam_per_unique_transcript$coverage)
                med_len_unique_tr <- median(length_per_unique_transcript$seqlengths)
                med_read_length <- median(bam_primary$readLength)

                # exporting 10% subsample primary data for statistical tests
                x <- nrow(bam_primary)
                subsample <- bam_primary[sample(1:nrow(bam_primary), round(x/10, digits=0)), ]

                subsample <- subsample %>%
                dplyr::mutate(sample=sample_name)
                if (!dir.exists(paste0(output,"csv/subsample"))){
                    dir.create(paste0(output,"csv/subsample"))}
                write.csv(subsample, file = paste0(output,"csv/subsample/",sample_name,".csv"), sep=",", quote=F, col.names = T, row.names=F)

                # Exporting the bam_primary dataframe (only primary alignment for each read)
                bam_primary <- bam_primary %>%
                dplyr::mutate(sample=sample_name)
                if (!dir.exists(paste0(output,"csv/primary_data"))){
                dir.create(paste0(output,"csv/primary_data"))}
                write.csv(bam_primary, file = paste0(output,"csv/primary_data/",sample_name,"_primary_data.csv"), sep=",", quote=F, col.names = T, row.names=F)

                # Exporting csv with each unique transcript and its associated length
                if (!dir.exists(paste0(output,"csv/per_unique_transcript"))){
                dir.create(paste0(output,"csv/per_unique_transcript"))}
                write.csv(bam_per_unique_transcript, file = paste0(output,"csv/per_unique_transcript/",sample_name,"_per_unique_transcript.csv"), sep=",", quote=F, col.names = T, row.names=F)

                
                # Creation of dataframes for proportions and numbers
                read_proportion_df <- bam_primary %>%
                dplyr::ungroup() %>%
                dplyr::select(novel_known, transcript_biotype) %>%
                dplyr::mutate(transcript_biotype = ifelse(transcript_biotype != "protein_coding" & transcript_biotype != "lncRNA", "others", as.character(transcript_biotype))) %>%
                as.data.frame()
                read_proportion_df[,"transcript_biotype"] <- as.factor(read_proportion_df[,"transcript_biotype"])

                read_proportion_df_plot <- read_proportion_df %>%
                dplyr::group_by(novel_known, transcript_biotype) %>%
                dplyr::summarise(count = n()) %>%
                dplyr::mutate(perc = (count / sum(count))) %>%
                dplyr::mutate(perc2 = round(count / sum(count) * 100, 1))

                read_nbr_df <- read_proportion_df %>%
                dplyr::select(novel_known) %>%
                dplyr::group_by(novel_known) %>%
                dplyr::summarise(count = n())

                # Ploting numbers
                nbr_plot_read <- ggplot(data=read_nbr_df, aes(x=novel_known,y=count)) +
                geom_bar(stat="identity", fill="#588157") +
                geom_text(aes(label=count), vjust=-0.3, size=13) +
                ggtitle("Number of read primary \naligned Novel vs Known") +
                xlab("Status of transcript") +
                ylab("Read count") +
                theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
                    axis.text = element_text(size = 30),
                    axis.title = element_text(size = 30))
                
                # Ploting proportions
                prop_plot_read <- ggplot(data=read_proportion_df_plot, aes(x=novel_known, y=perc2, label=perc2, fill=transcript_biotype)) +
                geom_bar(stat = "identity") +
                geom_text(position=position_stack(vjust = 0.5), size = 13) +
                xlab("Status of transcript") +
                ylab("Read proportion") +
                ggtitle("Proportion (%) of read primary \n aligned Novel vs Known") +
                theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 35),
                    legend.title = element_text(size = 25),
                    legend.text = element_text(size = 25),
                    axis.text = element_text(size = 30),
                    axis.title = element_text(size = 30))
                

                if (!dir.exists(paste0(output,"plots"))){
                    dir.create(paste0(output,"plots"))}

                if (!dir.exists(paste0(output,"plots/",sample_name))){
                    dir.create(paste0(output,"plots/",sample_name))}

                pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_alignment_read_count.pdf"), width=8, height=12)
                print(nbr_plot_read)
                dev.off()

                pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_alignment_read_proportions.pdf"), width=11, height=12)
                print(prop_plot_read)
                dev.off()

                rm(read_proportion_df)
                rm(read_proportion_df_plot)
                rm(read_nbr_df)

                message("Full data primary alignments proportion plots and bam_primary file exported")

            }


            # Creation of vectors with stats for each subgroups 
            cov_med_vec[[count]] <- median(bam_primary$coverage)
            read_len_med_vec[[count]] <- median(bam_primary$readLength)
            prop_full_tr_vec[[count]] <- sum(bam_primary$coverage > 0.95)/nrow(bam_primary)
            align_frac_med_vec[[count]] <- median(bam_primary$alignedFraction)
            accuracy_med_vec[[count]] <- median(bam_primary$accuracy)
            sec_align_mean_vec[[count]] <- mean(as.numeric(bam_primary$nbrSecondary))
            qscore_med_vec[[count]] <- median((bam_primary$mean_qscore_template))

            # Subgroups order for the vectors
            if(status=="full"){
            subgroups_order[[count]] <- "full"
            }
            else{subgroups_order[[count]] <- paste0(status,"_",biotype)}
            


            # Subgroups plots

            # Histogram of coverage
            
            plot1 <- ggplot(data=bam_primary, aes(x=coverage, fill=above)) +
            geom_histogram(bins = 180, show.legend = FALSE) +
            geom_vline(aes(xintercept=0.95), color="black", linetype="dashed", size=0.5) +
            theme_classic(base_size=16) +
            (if (status != "full"){
            ggtitle(paste0(status,"_",biotype))}
            else {ggtitle("full")}) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
            xlim(0.5,1) +
            xlab("Coverage Fraction") +
            ylab("Count") +
            scale_fill_manual(values = c("gray", "steelblue3"))
            plots1[[count]] <- plot1

            # Histogram of coverage vs length
            
            plot2 <- ggplot() +
            geom_hex(data=bam_primary, aes(x=seqlengths, y=coverage, fill = stat(log(count))), bins=100) +
            stat_smooth(data=bam_primary, aes(x=seqlengths, y=coverage), color="lavender", se=TRUE, size=0.5, level=0.95) +
            theme_classic(base_size=8) +
            (if (status != "full"){
            ggtitle(paste0(status,"_",biotype))}
            else {ggtitle("full")}) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            text = element_text(size = 16)) +
            xlim(0,15000) +
            ylim(0,1) +
            xlab("Transcript Length") +
            ylab("Coverage Fraction") +
            scale_fill_viridis_c()
            plots2[[count]] <- plot2


            # Secondary alignments bar chart
            
            plot3 <- ggplot(alignments) +
            geom_bar(stat='identity', aes(x=nbrSecondary, y=prop), fill = "steelblue3") +
            theme_classic(base_size=16) +
            (if (status != "full"){
            ggtitle(paste0(status,"_",biotype))}
            else {ggtitle("full")}) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
            xlab("Number of Secondary Alignments") +
            ylab("Proportion of Reads") +
            ylim(0,1) +
            scale_x_discrete(breaks = alignments$nbrSecondary, labels = alignments$nbrSecondary)
            plots3[[count]] <- plot3

            # Histogram unique transcript lengths
            
            plot4 <- ggplot(data=length_per_unique_transcript, aes(x=seqlengths)) +
            geom_histogram(bins = 180, show.legend = FALSE, fill="steelblue3") +
            theme_classic(base_size=16) +
            (if (status != "full"){
            ggtitle(paste0(status,"_",biotype))}
            else {ggtitle("full")}) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
            xlim(0,10000) +
            xlab("Transcript Length") +
            ylab("Count")
            plots4[[count]] <- plot4

            # Histogram of accuracy
            
            plot5 <- ggplot(data=bam_primary, aes(x=accuracy, y=..scaled..)) +
            geom_density(alpha = 0.4, show.legend = FALSE, fill="steelblue3") +
            theme_classic(base_size=16) +
            (if (status != "full"){
            ggtitle(paste0(status,"_",biotype))}
            else {ggtitle("full")}) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
            xlim(0.5,1) +
            xlab("Accuracy") +
            ylab("Density")
            plots5[[count]] <- plot5


            if (status=="full"){
                count <- count + 1
                message(paste0("Plots created for full data"))
                break
            } 

            count <- count + 1
            message(paste0("Plots created for subsample ",status,"_",biotype))
        }
    
    }


    # ---------------------------------------------------------------------------- #
    #             Exporting subgroups plots & general stats csv created            #
    # ---------------------------------------------------------------------------- #


    pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_coverage_fraction.pdf"), width=10, height=8)
    grid.arrange(grobs=plots1)
    dev.off()

    pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_density.pdf"), width=10, height=8)
    grid.arrange(grobs=plots2)
    dev.off()

    pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_sec_alns.pdf"), width=10, height=8)
    grid.arrange(grobs=plots3)
    dev.off() 

    pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_transcript_length_distribution.pdf"), width=10, height=8)
    grid.arrange(grobs=plots4)
    dev.off()

    pdf(paste0(output,"plots/",sample_name,"/",sample_name,"_accuracy.pdf"), width=10, height=8)
    grid.arrange(grobs=plots5)
    dev.off()

    message("Plots for subgroups exported")

    #Subgroups stats for general plots
    metric <- c("sample", #Sample
                "nbr_align", #Base number of alignements
                "nbr_reads", #Base number of unique reads
                "nbr_unique_tr", #Number of unique transcript
                "med_cov_unique_tr", #Unique transcripts coverage median 
                "med_len_unique_tr", #Unique transcripts length median
                "med_read_length", #Read length median
                "p_pcg_known", #Proportion of reads primary aligned to protein-coding known transcripts
                "p_pcg_novel", #Proportion of reads primary aligned to protein-coding novel transcripts
                "p_lnc_known", #Proportion of reads primary aligned to lncRNA known transcripts
                "p_lnc_novel", #Proportion of reads primary aligned to lncRNA novel transcripts
                "p_oth_known", #Proportion of reads primary aligned to other types known transcripts
                "p_oth_novel", #Proportion of reads primary aligned to other types novel transcripts
                paste0("cov_med_",subgroups_order[1]), #Transcript coverage median
                paste0("cov_med_",subgroups_order[2]),
                paste0("cov_med_",subgroups_order[3]),
                paste0("cov_med_",subgroups_order[4]),
                paste0("cov_med_",subgroups_order[5]),
                paste0("read_len_med_",subgroups_order[1]), #Read length median
                paste0("read_len_med_",subgroups_order[2]),
                paste0("read_len_med_",subgroups_order[3]),
                paste0("read_len_med_",subgroups_order[4]),
                paste0("read_len_med_",subgroups_order[5]),
                paste0("prop_full_tr_",subgroups_order[1]), #Prop of reads representing full length transcripts
                paste0("prop_full_tr_",subgroups_order[2]),
                paste0("prop_full_tr_",subgroups_order[3]),
                paste0("prop_full_tr_",subgroups_order[4]),
                paste0("prop_full_tr_",subgroups_order[5]),
                paste0("align_frac_med_",subgroups_order[1]), #Read aligned fraction median
                paste0("align_frac_med_",subgroups_order[2]),
                paste0("align_frac_med_",subgroups_order[3]),
                paste0("align_frac_med_",subgroups_order[4]),
                paste0("align_frac_med_",subgroups_order[5]),
                paste0("accuracy_med_",subgroups_order[1]), #Alignment accuracy median
                paste0("accuracy_med_",subgroups_order[2]),
                paste0("accuracy_med_",subgroups_order[3]),
                paste0("accuracy_med_",subgroups_order[4]),
                paste0("accuracy_med_",subgroups_order[5]),
                paste0("sec_align_mean_",subgroups_order[1]), #Secondary alignment mean
                paste0("sec_align_mean_",subgroups_order[2]),
                paste0("sec_align_mean_",subgroups_order[3]),
                paste0("sec_align_mean_",subgroups_order[4]),
                paste0("sec_align_mean_",subgroups_order[5]),
                paste0("qscore_med_",subgroups_order[1]), #Read mean qscore median
                paste0("qscore_med_",subgroups_order[2]),
                paste0("qscore_med_",subgroups_order[3]),
                paste0("qscore_med_",subgroups_order[4]),
                paste0("qscore_med_",subgroups_order[5])
                ) 
                
    outcome <- list(paste0(output,sample_name), nbr_align, nbr_read, nbr_unique_tr, med_cov_unique_tr, med_len_unique_tr, med_read_length,
                prop_read_known_pcg, prop_read_novel_pcg, prop_read_known_lnc, prop_read_novel_lnc, prop_read_known_oth, prop_read_novel_oth
                )

    for(i in 1:length(cov_med_vec)){
        outcome <- append(outcome, cov_med_vec[i])}
    for(i in 1:length(read_len_med_vec)){
        outcome <- append(outcome, read_len_med_vec[i])}
    for(i in 1:length(prop_full_tr_vec)){
        outcome <- append(outcome, prop_full_tr_vec[i])}
    for(i in 1:length(align_frac_med_vec)){
        outcome <- append(outcome, align_frac_med_vec[i])}
    for(i in 1:length(accuracy_med_vec)){
        outcome <- append(outcome, accuracy_med_vec[i])}
    for(i in 1:length(sec_align_mean_vec)){
        outcome <- append(outcome, sec_align_mean_vec[i])}
    for(i in 1:length(qscore_med_vec)){
        outcome <- append(outcome, qscore_med_vec[i])}

    # Metric as headers and outcome as single row values
    stats <- data.frame(matrix(ncol=length(metric), nrow=0))
    colnames(stats) <- metric
    stats[1,] <- outcome

        
    # Export overall alignments and reads proportions / stat file
    
    if (!dir.exists(paste0(output,"csv/sample_features"))){
            dir.create(paste0(output,"csv/sample_features"))}
    write.csv(stats, file = paste0(output,"csv/sample_features/",sample_name,".csv"), sep=",", quote=F, col.names = TRUE, row.names = FALSE) 

    message("Full data stats on bam_primary exported")

    }

    make_plots()

    message("Complete")
}

suppressWarnings(
  main())