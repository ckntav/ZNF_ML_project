# setwd("/Users/chris/Desktop/ZNF_ML_project")

library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(AnnotationDbi)
library(ef.utils)

#####
keepStdChr <- function(gr) {
  message("With all chromosomes, including contigs : ", length(gr), " regions")
  stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
  gr_StdChr <- keepSeqlevels(gr, stdChr[stdChr %in% seqlevels(gr)], pruning.mode = "coarse")
  message("Keeping standard chromosomes : ", length(gr_StdChr), " regions")
  message("\t--> ", length(gr) - length(gr_StdChr), " regions removed")
  return(gr_StdChr)
}

#
load_peak <- function(set_id) {
  message("# ", set_id)
  ENCODE_blacklist <- rtracklayer::import("input/ENCODE_blacklist_ENCFF419RSJ.bed")
  peaks_path <- file.path("output/peaks", paste0(set_id, ".bed"))
  peaks <- read_tsv(peaks_path, col_names = FALSE, col_types = "cddcdd") %>% dplyr::select(X1, X2, X3) %>% 
    set_names("seqnames", "start", "end") %>% GRanges
  message("Number of regions : ", length(peaks))
  if (sum(countOverlaps(peaks, ENCODE_blacklist)) != 0) {
    message("\t> Remove ", length(subsetByOverlaps(peaks, ENCODE_blacklist)), " blacklisted ENCODE regions")
    peaks <- subsetByOverlaps(peaks, ENCODE_blacklist, invert = TRUE)
    message("\tNumber of regions : ", length(peaks))
  }
  return(peaks)
}

raji_R1 <- load_peak("GSM3043267_Raji.R1") %>% keepStdChr()
raji_R2 <- load_peak("GSM3043268_Raji.R2") %>% keepStdChr()
u2os_R1 <- load_peak("GSM3043270_USOS.R1") %>% keepStdChr()
u2os_R2 <- load_peak("GSM3043271_USOS.R2") %>% keepStdChr()

gencode_v38 <- loadDb(file = "input/gencode/txdb.gencode38BasicChr.sqlite")

# annotate peaks
peakAnno_raji_R1 <- annotatePeak(raji_R1, TxDb=gencode_v38, annoDb="org.Hs.eg.db")
peakAnno_raji_R2 <- annotatePeak(raji_R2, TxDb=gencode_v38, annoDb="org.Hs.eg.db")
peakAnno_u2os_R1 <- annotatePeak(u2os_R1, TxDb=gencode_v38, annoDb="org.Hs.eg.db")
peakAnno_u2os_R2 <- annotatePeak(u2os_R2, TxDb=gencode_v38, annoDb="org.Hs.eg.db")

peakAllListAll = list("raji_R1" = peakAnno_raji_R1,
                      "raji_R2" = peakAnno_raji_R2,
                      "u2os_R1" = peakAnno_u2os_R1,
                      "u2os_R2" = peakAnno_u2os_R2)

plotAnnoBar(peakAllListAll)
plotDistToTSS(peakAllListAll, title = "Distribution of peaks relative to TSS")

write_tsv(peakAnno_raji_R1 %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_raji_R1.tsv")
write_tsv(peakAnno_raji_R2 %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_raji_R2.tsv")
write_tsv(peakAnno_u2os_R1 %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_u2os_R1.tsv")
write_tsv(peakAnno_u2os_R2 %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_u2os_R2.tsv")

# consensus peak raji
raji_list <- GRangesList("raji_R1" = raji_R1, "raji_R2" = raji_R2)
intersect_raji <- build_intersect(raji_list)
raji_twoR <- rowSums(intersect_raji$Matrix) >= 2
raji_consensus <- intersect_raji$Regions[raji_twoR]

peakAnno_raji_consensus <- annotatePeak(raji_consensus, TxDb=gencode_v38, annoDb="org.Hs.eg.db")
write_tsv(peakAnno_raji_consensus %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_raji_consensus.tsv")

# consensus peak u2os
u2os_list <- GRangesList("u2os_R1" = u2os_R1, "u2os_R2" = u2os_R2)
intersect_u2os <- build_intersect(u2os_list)
u2os_twoR <- rowSums(intersect_u2os$Matrix) >= 2
u2os_consensus <- intersect_u2os$Regions[u2os_twoR]

peakAnno_u2os_consensus <- annotatePeak(u2os_consensus, TxDb=gencode_v38, annoDb="org.Hs.eg.db")
write_tsv(peakAnno_u2os_consensus %>% as.data.frame,
          file = "output/peak_annotation/peakAnno_u2os_consensus.tsv")

#
peakAllConsensus <- list("raji_consensus" = peakAnno_raji_consensus,
                         "u2os_consensus" = peakAnno_u2os_consensus)
plotAnnoBar(peakAllConsensus)
plotDistToTSS(peakAllConsensus, title = "Distribution of peaks relative to TSS")
