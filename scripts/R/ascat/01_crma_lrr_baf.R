# ASCAT trial, step 1 - LRR and BAF from the CEL files with existing R tools.
#
# Nothing here is hand-written signal processing: aroma.affymetrix's CRMA v2
# (Bengtsson et al. 2009) does every step from CEL file to one total signal and one
# B-allele fraction per SNP:
#   allelic crosstalk calibration -> base position normalisation ->
#   probe summarisation (RmaCnPlm, the model aroma recommends for 10K-500K arrays) ->
#   PCR fragment length normalisation
# This script only turns its output into LRR:
#   LRR = log2(total / reference), 
#   with reference = median total of all arrays
# BAF is CRMA v2's allele B fraction as it comes.
#
# Run (repo root, ascat_env active):
#   Rscript scripts/R/ascat/01_crma_lrr_baf.R [cel_dir] [out_dir] [aroma_dir]
#   defaults: data/raw/cel  data/ascat/crma  data/ascat/aroma
# Writes (out_dir)   lrr.tsv.gz, baf.tsv.gz   probe | chr | pos | one column per array
#                    (same layout as data/processed/lrr.tsv.gz of step 06a; chr 23 = X, hg19)
#                    arrays.tsv               sample | gsm | x_level | n_X
# aroma_dir holds aroma's own folders (annotationData, rawData, probeData, plmData, ...).

suppressPackageStartupMessages({
  library(aroma.affymetrix)
  library(data.table)
})

CHIP <- "Mapping250K_Nsp"
DATASET <- "GSE44994"
ANNOTATION_URL <- "https://aroma-project.org/data/annotationData/chipTypes/Mapping250K_Nsp/"
ANNOTATION <- c("Mapping250K_Nsp.cdf",                      # which probe measures which SNP
                "Mapping250K_Nsp,na31,HB20101007.ugp",      # SNP positions (NetAffx na31 = hg19)
                "Mapping250K_Nsp,na31,HB20101007.ufl",      # PCR fragment lengths
                "Mapping250K_Nsp,HB20080710.acs")           # probe sequences

# aroma finds its files by folder name: annotationData/chipTypes/<chip>/ and rawData/<set>/<chip>/
prepare_aroma_folders <- function(cel_dir, aroma_dir) {
  annotation_dir <- file.path(aroma_dir, "annotationData", "chipTypes", CHIP)
  raw_dir <- file.path(aroma_dir, "rawData", DATASET, CHIP)
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in ANNOTATION) {
    target <- file.path(annotation_dir, f)
    if (file.exists(target)) next
    message("downloading ", f)
    zipped <- paste0(target, ".gz")
    download.file(paste0(ANNOTATION_URL, URLencode(basename(zipped))), zipped, mode = "wb", quiet = TRUE)
    R.utils::gunzip(zipped, target)
  }
  cels <- sort(list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE))
  stopifnot(length(cels) > 0)
  links <- file.path(raw_dir, basename(cels))
  missing <- !file.exists(links)
  if (any(missing)) file.symlink(normalizePath(cels)[missing], links[missing])
  basename(cels)
}

# X level of every array relative to all arrays: two groups, one copy and two copies
count_x_copies <- function(x_level) {
  groups <- kmeans(x_level, centers = range(x_level))
  ifelse(groups$cluster == which.max(groups$centers), 2L, 1L)
}

main <- function(args) {
  cel_dir   <- if (length(args) >= 1) args[1] else "data/raw/cel"
  out_dir   <- if (length(args) >= 2) args[2] else "data/ascat/crma"
  aroma_dir <- if (length(args) >= 3) args[3] else "data/ascat/aroma"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir)
  cel_files <- prepare_aroma_folders(cel_dir, aroma_dir)
  setwd(aroma_dir)

  # CEL files -> total signal and allele B fraction per SNP (CRMA v2)
  cdf <- AffymetrixCdfFile$byChipType(CHIP)
  cel_set <- AffymetrixCelSet$byName(DATASET, cdf = cdf)
  stopifnot(length(cel_set) == length(cel_files))
  message("CRMA v2 on ", length(cel_set), " arrays")
  crma <- doCRMAv2(cel_set, combineAlleles = FALSE, plm = "RmaCnPlm", verbose = -5) # computes baf (combine alleles set to false) and total, keeps them in crma$fracB, crma$total
  total <- extractMatrix(crma$total) # total signal
  frac_b <- extractMatrix(crma$fracB) # baf
  arrays <- getNames(crma$total) # "GSM1095396_SNP018ZXP2"
  stopifnot(identical(arrays, getNames(crma$fracB)))
  samples <- sub("^[^_]*_", "", arrays)
  colnames(total) <- colnames(frac_b) <- samples

  # SNPs on chromosomes 1-22 and X, in genome order
  positions <- getAromaUgpFile(crma$total)
  snps <- data.table(unit = seq_len(nbrOfUnits(cdf)), probe = getUnitNames(cdf),
                     chr = positions[, 1, drop = TRUE], pos = positions[, 2, drop = TRUE])
  snps <- snps[!is.na(chr) & chr %in% 1:23 & !is.na(pos) & startsWith(probe, "SNP_A-")]
  setorder(snps, chr, pos)
  total <- total[snps$unit, ]
  frac_b <- frac_b[snps$unit, ]
  total[!is.finite(total) | total <= 0] <- NA
  autosomal <- snps$chr <= 22

  # LRR against the median of all arrays, centred per array on its autosomes
  log_total <- log2(total)
  lrr <- log_total - matrixStats::rowMedians(log_total, na.rm = TRUE)
  lrr <- sweep(lrr, 2, matrixStats::colMedians(lrr[autosomal, ], na.rm = TRUE))

  # X: the reference must have two copies, so take it from the arrays with two X
  x_level <- matrixStats::colMedians(lrr[!autosomal, ], na.rm = TRUE)
  n_x <- count_x_copies(x_level) # kmeans on every array's X level relative to all arrays
  x_reference <- matrixStats::rowMedians(log_total[!autosomal, n_x == 2], na.rm = TRUE)
  centre <- matrixStats::colMedians(log_total[autosomal, ] -
                                    matrixStats::rowMedians(log_total[autosomal, ], na.rm = TRUE), na.rm = TRUE)
  lrr[!autosomal, ] <- sweep(log_total[!autosomal, ] - x_reference, 2, centre)

  frac_b[!is.finite(frac_b)] <- NA # if invalid then NaN, imagine a 0 diviion
  write_table <- function(m, name, digits) {
    out <- cbind(snps[, .(probe, chr, pos)], as.data.table(round(m, digits)))
    fwrite(out, file.path(out_dir, name), sep = "\t", na = "NA")
  }
  write_table(lrr, "lrr.tsv.gz", 4)
  write_table(pmin(pmax(frac_b, 0), 1), "baf.tsv.gz", 4) # ASCAT expects BAF within 0-1, so we clip it
  fwrite(data.table(sample = samples, gsm = sub("_.*$", "", arrays),
                    x_level = round(x_level, 4), n_X = n_x),
         file.path(out_dir, "arrays.tsv"), sep = "\t")
  message(nrow(snps), " SNPs x ", length(samples), " arrays -> ", out_dir,
          " (", sum(n_x == 1), " arrays with one X, ", sum(n_x == 2), " with two)")
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
