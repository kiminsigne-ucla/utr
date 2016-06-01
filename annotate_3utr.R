# generate a list of 3' UTR sequences
Sys.time()
library(dplyr)
library(ggplot2)

setwd("~/Documents/projects/utr")
# annotation previously trimmed to only include gene, transcript, and UTR entries for faster loading
annot <- read.table('gencode.v24.annotation.trimmed.gtf', header=F, sep='\t', comment.char='#', stringsAsFactors = F)
# annot <- read.table('test.gtf', header=F, comment.char='#', sep='\t', stringsAsFactors = F)
colnames(annot) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'info')

# add id column
annot <- mutate(annot, id = c(1:nrow(annot)))

# extract gene_id and transcript id from info column
fields <- strsplit(annot$info, ';')
# for each element in field, grab the 1st entry. Format: 'gene_id ENSG00000223972.5'. Then for each element in this list,
# split by ' ' and grab the second element, to grab only the IDs
gene_ids <- unlist( lapply( strsplit( unlist(lapply(fields, `[[`, 1)), ' '), `[[`, 2))
# same logic as gene_ids, but transcript ids is the second field in fields. Also, there is a leading whitespace so the format is
# " transcript_id ENST00000450305.2" and we can grab the third element in this case
transcript_ids <- unlist( lapply( strsplit( unlist(lapply(fields, `[[`, 2)), ' '), `[[`, 3))
annot <- mutate(annot, gene_id = gene_ids, transcript_id = transcript_ids)
# remove large vectors
rm(list=c('gene_ids', 'transcript_ids'))
# transcript id is only the second field when the type is transcript or UTR, otherwise it is gene type.
# If element doesn't begin with ENST, replace with NA
annot <- mutate(annot, transcript_id = ifelse(grepl('ENST', transcript_id), transcript_id, NA))

# how many transcripts per gene?
transcript_counts <- annot %>% 
                    group_by(gene_id) %>%
                    summarise(count = n_distinct(transcript_id, na.rm=T))
table(transcript_counts$count)

# UTRs per transcript
utr_counts <- annot %>%
            group_by(transcript_id) %>%
            summarise(utr_count = sum(type == 'UTR'))
table(utr_counts$utr_count)

is_3UTR <- function(transcript_start, transcript_end, strand, utr_start, utr_end) {
    # Define as 3' UTR if it falls in the second half of the transcript
    midpoint <- transcript_start + ( (transcript_end - transcript_start)/2 )
    if(strand == '+'){
        if(utr_start > midpoint){ return(TRUE) }
        else{ return(FALSE) }
    }
    if(strand == '-'){
        if(utr_start < midpoint){ return(TRUE) }
        else{ return(FALSE) }
    }
}

# # deal with transcripts that have 0 UTRs separately
# no_utr_ids <- filter(utr_counts, utr_count == 0) %>% select(transcript_id)
# no_utrs <- filter(annot, transcript_id %in% no_utr_ids$transcript_id)
# no_utrs <- mutate(no_utrs, is_3utr = rep(NA, nrow(no_utrs)))

transcript_with_utr_ids <- filter(utr_counts, utr_count > 0) %>% select(transcript_id)
transcript_with_utr <- filter(annot, transcript_id %in% transcript_with_utr_ids$transcript_id)
# initialize empty column first, faster to access by index and change than it is to build and append new vector
transcript_with_utr <- mutate(transcript_with_utr, is_3utr = rep(NA, nrow(transcript_with_utr)))

Sys.time()

# this part takes a while but only needs to be run once to generate the 3' UTR annotation
for(transcript in unique(transcript_with_utr$transcript_id)){
    info <- filter(transcript_with_utr, transcript_id == transcript)
    transcript_info <- filter(info, type == 'transcript')
    transcript_start <- transcript_info$start[1]
    transcript_end <- transcript_info$end[1]
    strand <- transcript_info$strand[1]
    
    utrs <- filter(info, type == 'UTR')
    for(i in 1:nrow(utrs)){
        utr_start <- utrs$start[i]
        utr_end <- utrs$end[i]
        is_3utr <- is_3UTR(transcript_start = transcript_start, transcript_end = transcript_end, strand = strand,
                           utr_start = utr_start, utr_end = utr_end)
        index <- which(transcript_with_utr$id == utrs$id[i])
        # subset original data frame and change is_3utr value
        transcript_with_utr$is_3utr[index] = is_3utr
    }
}
Sys.time()
# write those that are 3' UTR to file
is_3utrs <- filter(transcript_with_utr, is_3utr == TRUE) %>% select(chr, start, end, strand, gene_id, transcript_id)
write.table(is_3utrs, file='annotated_3UTR.tsv', col.names=T, row.names=F, quote=F, sep='\t')

