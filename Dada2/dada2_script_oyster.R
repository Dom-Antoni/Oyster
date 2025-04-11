if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
.libPaths()
BiocManager::install("dada2", dependencies = T)
cat('.libPaths("~/Rlibs")', file = "~/.Rprofile", append = TRUE)

install.packages(c("BiocManager", "devtools", dependencies = T))
BiocManager::version()
BiocManager::install("ShortRead")

install.packages("dplyr")
install.packages("DECIPHER")
install.packages("tidyverse")

BiocManager::install("DECIPHER")


library("dada2")
library("Rcpp")
library("ggplot2")
library("dplyr")
library("ShortRead")
library("tidyverse")

setwd("/home/dantoni/Desktop/Austern/ASV_Fasta/Wormhole m92K6z")
samples <- scan("samples", what = "Character")

forward_reads <- paste0(samples, "_sub_R1.fastq.gz")
reverse_reads <- paste0(samples, "_sub_R2.fastq.gz")

#----
#This code chung just analyses what reads are most frequently occuring in the raw sequences 
#of the gill samples. This codes does not need to be run to reproduce the findings but I dont want to
# delte it since it might be interessting for others. 
get_most_frequent_read <- function(file) {
  fq <- readFastq(file)
  seqs <- as.character(sread(fq))
  seq_table <- table(seqs)
  most_frequent <- names(which.max(seq_table))
  freq <- max(seq_table)
  sample_name <- basename(file)
  
  return(data.frame(
    Sample = sample_name,
    MostFrequentRead = most_frequent,
    Frequency = freq,
    TotalReads = length(seqs)
  ))
}

results_f <- map_dfr(forward_reads, get_most_frequent_read)
results_r <- map_dfr(reverse_reads, get_most_frequent_read)

results_f <- results_f %>%
  mutate(Sample = as.numeric(str_extract(Sample, "^[^\\-]+")))%>%
  arrange(Sample)

results_r <- results_r %>%
  mutate(Sample = as.numeric(str_extract(Sample, "^[^\\-]+")))%>%
  arrange(Sample)

Gill_read <- results_f$MostFrequentRead[1]
Gill_read_r <- results_r$MostFrequentRead[1]

results_f$MostFrequentRead == Gill_read
results_r$MostFrequentRead == Gill_read_r

median(results_f$TotalReads[1:45])
median(results_f$TotalReads[45:87])


ggplot(results_r, aes(x= Sample, y= TotalReads))+
  geom_col()+
  theme_minimal()
  
###----

filtered_forward_reads <- paste0(samples,"_filtered_forward.fastq.gz")
filtered_reverse_reads <- paste0(samples,"_filtered_reverse.fastq.gz")

p_quality_f_reads_unfiltered <- plotQualityProfile(forward_reads)
p_quality_r_reads_unfiltered <- plotQualityProfile(reverse_reads)

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, compress = T,
                              reverse_reads, filtered_reverse_reads, maxEE = c(2,2),
                              rm.phix = T, minLen = 170, truncLen= c(250,200))

filtered_out_per <- filtered_out %>%
  as.data.frame() %>%
  mutate(percent = reads.out *100 / reads.in )

results_f_filtered <- map_dfr(filtered_forward_reads, get_most_frequent_read)
results_r_filtered <- map_dfr(filtered_reverse_reads, get_most_frequent_read)

results_f_filtered <- results_f_filtered %>%
  mutate(Sample = as.numeric(str_extract(Sample, "^[^\\-]+")))%>%
  arrange(Sample)%>%
  mutate(frac = Frequency / TotalReads * 100)

median(results_f_filtered$frac[1:45])
median(results_f_filtered$frac[45:87])


results_r_filtered <- results_r_filtered %>%
  mutate(Sample = as.numeric(str_extract(Sample, "^[^\\-]+")))%>%
  arrange(Sample)



filtered_out_2 <- filterAndTrim(forward_reads, filtered_forward_reads, compress = T,
                              reverse_reads, filtered_reverse_reads, maxEE = c(2,2),
                              rm.phix = T, minLen = 170, truncLen= c(225,180))

p_quality_f_reads_filtered <- plotQualityProfile(filtered_forward_reads) 
p_quality_r_reads_filtered <- plotQualityProfile(filtered_reverse_reads) 

p_quality_f_reads_filtered2 <- plotQualityProfile(filtered_forward_reads) 
p_quality_r_reads_filtered2 <- plotQualityProfile(filtered_reverse_reads) 

# I chose the more stringend filtering step with 225 and 180 as TruncLen 
# It just has better quality 

err_forward_reads <- learnErrors(filtered_forward_reads)
err_reverse_reads <- learnErrors(filtered_reverse_reads)

plotErrors(err_forward_reads, nominalQ=T)
plotErrors(err_reverse_reads, nominalQ=T)

derep_forward <- derepFastq(filtered_forward_reads, verbose = T)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = T)
names(derep_reverse) <- samples

dada_forward <- dada(derep_forward, err= err_forward_reads, pool=T)
dada_reverse <- dada(derep_reverse, err= err_reverse_reads, pool=T)

merged_amplicons <- mergePairs(dada_forward, derep_forward, 
                               dada_reverse, derep_reverse, trimOverhang = T, minOverlap = 12)


merged_amplicons$`10-K-T1-Oli1-B_S12__filtered_forward.fastq.gz`

class(merged_amplicons)
length(merged_amplicons)
names(merged_amplicons)

seqtab <- makeSequenceTable(merged_amplicons)

dim(seqtab)

seqtab_nochim <- removeBimeraDenovo(seqtab, verbose = T) 

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out_2[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab_nochim),
                          final_perc_reads_retained=round(rowSums(seqtab_nochim)/filtered_out_2[,1]*100, 1))

write.table(summary_tab, "read-count-tracking_oyster.tsv", quote=FALSE, sep="\t", col.names=NA)

load("SILVA_SSU_r138_2019.RData")
library("DECIPHER")
packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab_nochim))

tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand= "both", processors = NULL)

load("tax-info.RData")

asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")
for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))

write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab_nochim)
row.names(asv_tab) <- sub(">","",asv_headers)

write.table(asv_tab,"ASVs_count_table.tsv", sep=",", quote = F, col.names = NA)


ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxoonomy.tsv", sep = ",", quote=F, col.names=NA)
            