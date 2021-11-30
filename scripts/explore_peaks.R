# setwd("/Users/chris/Desktop/ZNF_ML_project")

library(tidyverse)
library(GenomicRanges)

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




#####
load_cofactor_peaks <- function(cofactors = c("MED1", "BRD4", "CDK9", "NIPBL", "SMC1A"), type) {
  ENCODE_blacklist <- rtracklayer::import("input/ENCODE_blacklist_ENCFF419RSJ.bed")
  peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"
  # artefact <- GRanges(seqnames = "chr3", IRanges(start = 93470340, end = 93470809))
  cofactors_peaks <- GRangesList()
  name_cofactors_peaks <- c()
  for (cofactor in cofactors) {
    for (condition in c("CTRL", "DEX")) {
      message("####\t", cofactor, " | ", condition)
      basename <- paste("A549", condition,cofactor, "rep1", type, sep = "_")
      peaks_path <- file.path(peaks_dir, basename, paste0(basename, "_peaks.", type, "Peak.bed"))
      # message(peaks_path)
      peaks <- rtracklayer::import(peaks_path)
      message("Number of regions : ", length(peaks))
      if (sum(countOverlaps(peaks, ENCODE_blacklist)) != 0) {
        message("\t> Remove ", length(subsetByOverlaps(peaks, ENCODE_blacklist)), " blacklisted ENCODE regions")
        peaks <- subsetByOverlaps(peaks, ENCODE_blacklist, invert = TRUE)
        message("\tNumber of regions : ", length(peaks))
      }
      cofactors_peaks <- append(cofactors_peaks, GRangesList(peaks))
      name_cofactors_peaks <- c(name_cofactors_peaks, paste0(cofactor, "_", condition))
    }
  }
  names(cofactors_peaks) <- name_cofactors_peaks
  message("#####################################")
  message("Available set of regions: ")
  print(names(cofactors_peaks))
  return(cofactors_peaks)
}


