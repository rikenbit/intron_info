
################################
# Load libraries
library(tidyverse)
library(rtracklayer)

# Make output directory
if(! dir.exists(path_outdir)){
	dir.create(path_outdir, recursive=TRUE)
}

# Read GTF file
options(timeout = 300)
z <- import(path_gtf)

# Select `exon` entries in GTF
z2 <- z[z$type == "exon"]

# Make tibble
dfgtf <- tibble(
	transcript_id = z2$transcript_id, 
	gene_type = z2$gene_type,
	transcript_type = z2$transcript_type,
	exon_length = width(z2),
	exon_start = start(z2)
	)

dfgtf %>% 
	mutate(transcript_support_level = z2$transcript_support_level) -> dfgtf

# 
dfgtf %>%
	arrange(transcript_id, exon_start) %>%
	group_by(transcript_id, gene_type, transcript_type, transcript_support_level) %>%
	summarise(
		exon_num = n(),
		exon_length_list = list(exon_length),
		exon_start_list = list(exon_start)
		) -> dfgtf

# Calculate intron length
calc_intron_length <- function(exon_length_list, exon_start_list){
	sapply(1:(length(exon_length_list)-1), function(i){
		as.integer(exon_start_list[i+1])-as.integer(exon_start_list[i])-as.integer(exon_length_list[i])
	})
}

dfgtf %>%
	mutate(
		intron_length_list = map2(exon_length_list, exon_start_list, calc_intron_length)
	) -> dfgtf


# Stat by biotype
dfgtf %>%
	filter(exon_num > 1) %>%
	group_by(gene_type) %>%
	summarise(
		intron_len_min=min(unlist(.data[["intron_length_list"]])),
		intron_len_q1=quantile(unlist(.data[["intron_length_list"]]), probs=0.25),
		intron_len_median=median(unlist(.data[["intron_length_list"]])),
		intron_len_q3=quantile(unlist(.data[["intron_length_list"]]), probs=0.75),
		intron_len_max=max(unlist(.data[["intron_length_list"]])),
		n_transcript = n()
	) -> dfgtfsummary

dfgtfsummary %>% arrange(-n_transcript) -> dfgtfsummary

# Stat by biotype
dfgtf %>%
	filter(exon_num > 1) %>%
	group_by(gene_type, transcript_type) %>%
	summarise(
		intron_len_min=min(unlist(.data[["intron_length_list"]])),
		intron_len_q1=quantile(unlist(.data[["intron_length_list"]]), probs=0.25),
		intron_len_median=median(unlist(.data[["intron_length_list"]])),
		intron_len_q3=quantile(unlist(.data[["intron_length_list"]]), probs=0.75),
		intron_len_max=max(unlist(.data[["intron_length_list"]])),
		n_transcript = n()
	) -> dfgtfsummary2

dfgtfsummary2 %>% arrange(-n_transcript) -> dfgtfsummary2

# Stat exon number
dfgtf %>%
	group_by(gene_type) %>%
	summarise(
		exon_num_min=min(exon_num),
		exon_num_q1=quantile(exon_num, probs=0.25),
		exon_num_median=median(exon_num),
		exon_num_q3=quantile(exon_num, probs=0.75),
		exon_num_max=max(exon_num),
		n_transcript = n()
	) %>% arrange(-n_transcript) -> dfsummary_exon_num
dfgtf %>%
	group_by(gene_type, transcript_type) %>%
	summarise(
		exon_num_min=min(exon_num),
		exon_num_q1=quantile(exon_num, probs=0.25),
		exon_num_median=median(exon_num),
		exon_num_q3=quantile(exon_num, probs=0.75),
		exon_num_max=max(exon_num),
		n_transcript = n()
	) %>% arrange(-n_transcript) -> dfsummary_exon_num2

#

# Min intron length per transcript
dfgtf %>%
	filter(exon_num > 1) %>%
	mutate(min_intron_length_per_transcript = map_dbl(intron_length_list, min)) %>%
	select(transcript_id, gene_type, transcript_type, transcript_support_level, min_intron_length_per_transcript) -> dfgtf_min_intron_length_per_transcript

# Count transcripts with short introns
dfgtf_min_intron_length_per_transcript %>%
	group_by(gene_type, transcript_type, transcript_support_level) %>%
	summarise(
		n_min10 = sum(min_intron_length_per_transcript <= 10),
		n_min20 = sum(min_intron_length_per_transcript <= 20),
		n_min30 = sum(min_intron_length_per_transcript <= 30),
		n_min40 = sum(min_intron_length_per_transcript <= 40),
		n_min50 = sum(min_intron_length_per_transcript <= 50),
		n_transcript = n()
		) %>% 
	arrange(-n_transcript) -> dfgtf_min_intron_length_per_transcript_count

# Save
write_tsv(dfgtfsummary,
	file.path(path_outdir, "intron_length_summary.tsv"))
write_tsv(dfgtfsummary2,
	file.path(path_outdir, "intron_length_summary_transcript_type.tsv"))

write_tsv(dfsummary_exon_num,
	file.path(path_outdir, "exon_number_summary.tsv"))
write_tsv(dfsummary_exon_num2,
	file.path(path_outdir, "exon_number_summary_transcript_type.tsv"))

write_tsv(dfgtf_min_intron_length_per_transcript,
	file.path(path_outdir, "min_intron_length_per_transcript.tsv"))

write_tsv(dfgtf_min_intron_length_per_transcript_count,
	file.path(path_outdir, "count_min_intron_length_per_transcript.tsv"))


save(dfgtf, dfgtfsummary, dfgtf_min_intron_length_per_transcript,
	dfgtfsummary2, dfsummary_exon_num, dfsummary_exon_num2,
	dfgtf_min_intron_length_per_transcript_count, 
	file=file.path(path_outdir, "intron_length_master_table.Rdata"))

# SessionInfo()
sessionInfo()


