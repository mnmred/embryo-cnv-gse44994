# Step 06b - cut every chromosome of every array into segments with copynumber::aspcf.
#
# ASPCF (allele-specific piecewise constant fitting, Nilsen et al. 2012) segments LRR and
# BAF together: a breakpoint can come from a change in the amount of DNA (LRR) or in the
# balance of the two alleles (BAF). It was made for mixtures of cells (tumour samples),
# which also suits embryo biopsies where only some cells may carry a change.
# Settings: copynumber's defaults (gamma = 40; LRR winsorized first).
#
# ASPCF only cuts the genome into pieces and reports each piece's mean LRR and BAF.
# Step 06c (scripts/python/06c_call_copy_number.py) gives every piece a copy number.
#
# Run (repo root, bioenv active):
#   Rscript scripts/R/06b_segment.R [out_dir] [cores]      defaults: data/processed 4
#   about 4 min with 6 cores
# Reads   (out_dir)  lrr.tsv.gz, baf.tsv.gz   (step 06a)
# Writes  (out_dir)  segments_raw_ASPCF.tsv
#                    sample | chr | start | end (bp of the first / last SNP) | n_snps |
#                    mean (mean LRR of the piece) | mean_baf
# File names must match snparray/paths.py (lrr_tsv, baf_tsv, segments_raw).


suppressPackageStartupMessages({
  library(copynumber)
  library(data.table)
  library(parallel)
})

# chrom | pos | one LRR (or BAF) column named after the array: the input copynumber expects
one_array <- function(table, s) {
  cbind(data.frame(chrom = table$chr, pos = table$pos), setNames(data.frame(table[[s]]), s))
}

# ASPCF on one array: pull extreme LRR values in (winsorize), then segment LRR + BAF
segment_aspcf <- function(lrr, baf, s, gamma = 40) {
  w <- winsorize(one_array(lrr, s), assembly = "hg19", verbose = FALSE)
  out <- aspcf(w, one_array(baf, s), gamma = gamma, assembly = "hg19", verbose = FALSE)
  data.table(sample = s, chr = out$chrom, start = out$start.pos, end = out$end.pos,
             n_snps = out$n.probes, mean = out$logR.mean, mean_baf = out$BAF.mean)
}

main <- function(args) {
  out_dir <- if (length(args) >= 1) args[1] else "data/processed"
  cores   <- if (length(args) >= 2) as.integer(args[2]) else 4

  # get the LRR and BAF tables, one column per array, one row per SNP
  lrr <- fread(file.path(out_dir, "lrr.tsv.gz"))
  baf <- fread(file.path(out_dir, "baf.tsv.gz"))
  samples <- setdiff(names(lrr), c("probe", "chr", "pos"))
  message("ASPCF on ", length(samples), " arrays x ", nrow(lrr), " SNPs, ", cores, " cores")

  t0 <- Sys.time()
  res <- mclapply(samples, function(s) segment_aspcf(lrr, baf, s), mc.cores = cores)
  # errors can happen if a sample has no SNPs on a chromosome, or if all LRR values are extreme and get winsorized to the same value
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) stop("ASPCF failed for ", paste(samples[failed], collapse = ", "), ": ", res[[which(failed)[1]]])
  seg <- rbindlist(res)

  target <- file.path(out_dir, "segments_raw_ASPCF.tsv")
  fwrite(seg, target, sep = "\t")
  message(nrow(seg), " segments -> ", target, " (",
          format(round(difftime(Sys.time(), t0, units = "secs"))), ")")
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
