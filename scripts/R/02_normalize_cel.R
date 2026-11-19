# Step 02 - read the CEL files and normalise probe intensities.
#
# Writes (into out_dir)
#   pm_log2.f32         normalised log2 PM intensities, float32, one block per array
#   pm_probes.tsv       per probe: fid, fsetid (SNP id), allele (0 = A, 1 = B), strand
#   snp_annotation.tsv  per SNP: id, name, chromosome, hg19 position, cytoband,
#                       PCR fragment length, probe GC fraction
#   cel_order.txt       array order of the blocks in pm_log2.f32

suppressPackageStartupMessages({
  library(oligo)
  library(pd.mapping250k.nsp)
  library(DBI)
})

# CEL files -> matrix of raw PerfectMatch intensities (probes x cels), rownames = probe fid
read_pm_intensities <- function(cel_dir) {
  # read all CEL files, sort by name (GSM id) to ensure reproducible order, and extract the PM matrix
  cels <- sort(list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE))
  stopifnot(length(cels) > 0) # checks files were found (preventive)
  message("reading ", length(cels), " CEL files from ", cel_dir)
  raw <- read.celfiles(cels, verbose = FALSE) # oligo
  pm_matrix <- pm(raw) # keep only the PerfectMatch probes (ignore the Mismatch probes)
  attr(pm_matrix, "cel_files") <- basename(cels)
  pm_matrix
}

# Quantile normalisation: every cel gets exactly the reference distribution
# a probe keeps its rank within its cels (ties share the mean reference value)
quantile_normalize_to_reference <- function(pm_matrix) {
  # create variables stored in the pd.mapping250k.nsp package (the reference distribution is the HapMap target)
  load(system.file("extdata", "pd.mapping250k.nspRef.rda", package = "pd.mapping250k.nsp"))
  target <- sort(reference)
  n <- nrow(pm_matrix)
  stopifnot(length(target) == n) # safety check: the reference distribution has one reference per probe
  for (j in seq_len(ncol(pm_matrix))) {
    r <- rank(pm_matrix[, j], ties.method = "average")
    pm_matrix[, j] <- approx(seq_len(n), target, xout = r)$y # use approx to interpolate the reference value for each rank and take care of the ties
  }
  pm_matrix
}

# Which SNP / allele / strand each probe row measures
probe_table <- function(fid) {
  con <- db(pd.mapping250k.nsp)
  p <- dbGetQuery(con, "select fid, fsetid, allele, strand from pmfeature")
  p <- p[match(fid, p$fid), ]
  stopifnot(identical(p$fid, fid)) # safety check: probes re ordering worked
  p
}

# SNP positions (hg19), cytobands, PCR fragment length and probe GC content
snp_annotation <- function() {
  con <- db(pd.mapping250k.nsp)
  ann <- dbGetQuery(con, "select fsetid, man_fsetid, chrom, physical_pos, cytoband from featureSet")
  fl  <- dbGetQuery(con, "select fsetid, max(length) frag_len from fragmentLength group by fsetid")
  sq  <- dbGetQuery(con, "select p.fsetid, s.seq from pmfeature p join sequence s on p.fid = s.fid")
  sq$gc <- nchar(gsub("[AT]", "", sq$seq)) / nchar(sq$seq)
  gc  <- aggregate(gc ~ fsetid, sq, mean)
  merge(merge(ann, fl, all.x = TRUE), gc, all.x = TRUE) # left join to keep all SNPs, even if some have no fragment length or GC content
}

write_probe_level <- function(pm_norm, probes, cel_files, annotation, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # log2 so brightness differences are additive (log2(A/B) = log2(A) - log2(B))
  writeBin(as.numeric(log2(pm_norm)), file.path(out_dir, "pm_log2.f32"), size = 4)  # column-by-column float32, one block per cel
  write.table(probes, file.path(out_dir, "pm_probes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(annotation, file.path(out_dir, "snp_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(cel_files, file.path(out_dir, "cel_order.txt"))
}

# Standard deviation of every column (= every array), computed one column at a time
# so the whole matrix is never copied. `transform` is applied to the column first.
column_sd <- function(m, transform = identity) {
  vapply(seq_len(ncol(m)), function(j) sd(transform(m[, j])), numeric(1))
}

main <- function(args) {
  cel_dir <- if (length(args) >= 1) args[1] else "data/raw/cel"
  out_dir <- if (length(args) >= 2) args[2] else "data/processed/probe_level"
  pm_raw <- read_pm_intensities(cel_dir)
  message("read ", ncol(pm_raw), " arrays x ", nrow(pm_raw), " PM probes")
  cel_files <- attr(pm_raw, "cel_files")
  fid <- as.integer(rownames(pm_raw))
  sd_raw <- column_sd(pm_raw)          # before normalisation (brightness units)
  sd_raw_log2 <- column_sd(pm_raw, log2)    # before normalisation, on the log2 scale
  pm_norm <- quantile_normalize_to_reference(pm_raw)
  sd_norm     <- column_sd(pm_norm)         # after quantile normalisation
  sd_norm_log2 <- column_sd(pm_norm, log2)  # after normalisation + log2 (what gets saved)
  sd_table <- data.frame(array = cel_files,
                         raw = sd_raw, raw_log2 = sd_raw_log2,
                         quantile = sd_norm, quantile_log2 = sd_norm_log2)
  print(sd_table, digits = 4)
  rm(pm_raw); invisible(gc())
  write_probe_level(pm_norm, probe_table(fid), cel_files, snp_annotation(), out_dir)
  message("written to ", out_dir)
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
