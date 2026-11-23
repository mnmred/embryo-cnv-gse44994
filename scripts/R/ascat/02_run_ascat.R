# ASCAT trial, step 2 - copy number of every segment with ASCAT (Van Loo et al. 2010).
#
# ASCAT does both jobs of step 06 of the main pipeline in one package:
#   ascat.predictGermlineGenotypes  which SNPs are homozygous (no matched normal needed;
#                                   ASCAT ships settings for this array: "Affy250k_nsp")
#   ascat.aspcf                     cuts LRR and BAF into segments (ASPCF, penalty 70)
#   ascat.runAscat                  gives every segment an integer copy number per allele
#
# ASCAT was made for tumours, so it has two parameters an embryo biopsy does not need,
# and one that depends on the array. They are tried in three settings:
#   default      gamma 0.55, purity and ploidy fitted by ASCAT       (ASCAT as it comes)
#   fixed        gamma 0.55, purity 1 and ploidy 2 given             (an embryo is not a mixture)
#   calibrated   gamma from the X of the one-X arrays, purity 1, ploidy 2
# gamma is the LRR shift of one lost copy (1 in theory; ASCAT's default for arrays is 0.55).
#
# Run (repo root, ascat_env active):
#   Rscript scripts/R/ascat/02_run_ascat.R <in_dir> <out_dir> [cores]
#   pipeline LRR/BAF:  Rscript scripts/R/ascat/02_run_ascat.R data/processed   data/ascat/from_pipeline
#   CRMA v2 LRR/BAF:   Rscript scripts/R/ascat/02_run_ascat.R data/ascat/crma  data/ascat/from_crma
# Reads   (in_dir)   lrr.tsv.gz, baf.tsv.gz   probe | chr | pos | one column per array (chr 23 = X)
#                    arrays.tsv               sample, n_X (number of X copies)
# Writes  (out_dir)  segments_ASCAT_<setting>.tsv   sample | chr | start | end | n_major | n_minor | copies
#                    fit_ASCAT.tsv                  per array and setting: purity, ploidy, goodness of fit, status
#                    gamma.tsv                      the calibrated gamma
#                    work/                          ASCAT's input files, segmentation and plots

suppressPackageStartupMessages({
  library(ASCAT)
  library(data.table)
  library(parallel)
})

PLATFORM <- "Affy250k_nsp" # name for this chip in ASCAT
DEFAULT_GAMMA <- 0.55
SEED <- 1

# ASCAT reads text tables: row names = SNP, then chrs, pos, one column per array
write_ascat_input <- function(table, samples, target) {
  out <- data.frame(chrs = ifelse(table$chr == 23, "X", as.character(table$chr)),  # ASCAT wants "X" not 23
                    pos = as.integer(table$pos), table[, samples, with = FALSE],
                    row.names = table$probe, check.names = FALSE) # don't change the sample names
  write.table(out, target, sep = "\t", quote = FALSE, col.names = NA)
}

# LRR shift of one lost copy: the X of arrays with one X against the X of arrays with two (male vs female)
gamma_from_x <- function(lrr, samples, n_x) {
  x_level <- vapply(samples, function(s) median(lrr[chr == 23][[s]], na.rm = TRUE), numeric(1)) #median lrr over all X SNPs
  median(x_level[n_x == 2]) - median(x_level[n_x == 1]) # median over arrays with two X minus median over arrays with one X
}

# one array at a time, so a failure of one array does not stop the others
run_one <- function(segmented, i, gamma, fixed, plot_dir) {
  one <- segmented
  keep <- function(m) m[, i, drop = FALSE]
  for (name in c("Tumor_LogR", "Tumor_BAF", "Tumor_LogR_segmented")) one[[name]] <- keep(one[[name]])
  one$Tumor_BAF_segmented <- one$Tumor_BAF_segmented[i] #a list, not matrix
  one$samples <- one$samples[i]; one$gender <- one$gender[i]
  s <- one$samples
  out <- tryCatch(
    suppressWarnings(ascat.runAscat(one, gamma = gamma, write_segments = FALSE, img.dir = plot_dir,
                                    rho_manual = if (fixed) 1 else NA, psi_manual = if (fixed) 2 else NA)),
    error = function(e) NULL)
  ok <- !is.null(out) && !is.null(out$segments) && nrow(out$segments) > 0
  fit <- data.table(sample = s,
                    purity = if (ok) unname(out$purity[1]) else NA_real_,
                    ploidy = if (ok) unname(out$ploidy[1]) else NA_real_,
                    goodness_of_fit = if (ok) unname(out$goodnessOfFit[1]) else NA_real_,
                    status = if (!ok) "no solution" else if (length(out$nonaberrantarrays) > 0) "no aberration" else "fitted")
  segments <- if (ok) as.data.table(out$segments)[, .(
    sample, chr = ifelse(chr == "X", 23L, suppressWarnings(as.integer(chr))),
    start = startpos, end = endpos, n_major = nMajor, n_minor = nMinor, copies = nMajor + nMinor)] else NULL
  list(fit = fit, segments = segments)
}

main <- function(args) {
  if (length(args) < 2) stop("usage: 02_run_ascat.R <in_dir> <out_dir> [cores]")
  in_dir <- args[1]; out_dir <- args[2]
  cores <- if (length(args) >= 3) as.integer(args[3]) else 4
  work <- file.path(out_dir, "work")
  dir.create(file.path(work, "plots"), recursive = TRUE, showWarnings = FALSE)

  # inputs: LRR and BAF in the same SNP order, plus the number of X copies of every array
  lrr <- fread(file.path(in_dir, "lrr.tsv.gz"))
  baf <- fread(file.path(in_dir, "baf.tsv.gz"))
  stopifnot(identical(lrr$probe, baf$probe))
  arrays <- fread(file.path(in_dir, "arrays.tsv"))
  samples <- setdiff(names(lrr), c("probe", "chr", "pos"))
  n_x <- arrays$n_X[match(samples, arrays$sample)]
  stopifnot(!anyNA(n_x))

  # gamma for the "calibrated" setting, written down next to the default for the record
  gamma_x <- gamma_from_x(lrr, samples, n_x)
  fwrite(data.table(setting = c("default", "fixed", "calibrated"),
                    gamma = c(DEFAULT_GAMMA, DEFAULT_GAMMA, round(gamma_x, 4))),
         file.path(out_dir, "gamma.tsv"), sep = "\t")
  message(length(samples), " arrays x ", nrow(lrr), " SNPs; gamma from X: ", round(gamma_x, 3))

  write_ascat_input(lrr, samples, file.path(work, "Tumor_LogR.txt"))
  write_ascat_input(baf, samples, file.path(work, "Tumor_BAF.txt"))
  # gender tells ASCAT whether X is expected to have 1 copy (XY) or 2 (XX)
  loaded <- ascat.loadData(file.path(work, "Tumor_LogR.txt"), file.path(work, "Tumor_BAF.txt"),
                           chrs = c(1:22, "X"), gender = ifelse(n_x == 2, "XX", "XY"),
                           genomeVersion = "hg19")

  # homozygous SNPs, then segmentation: the slow part, done once and kept
  segmented_file <- file.path(work, "segmented.rds")
  if (file.exists(segmented_file)) {
    segmented <- readRDS(segmented_file)
  } else {
    t0 <- Sys.time()
    # Without a matched normal (e.g. the parent's blood) ASCAT cannot know each SNP's genotype.
    # It guesses from BAF: SNPs near 0 or 1 are homozygous and carry no allele balance, so only
    # the remaining (heterozygous) SNPs are used for BAF. Noisy BAF makes this guess harder.
    germline <- ascat.predictGermlineGenotypes(loaded, platform = PLATFORM,
                                               img.dir = file.path(work, "plots"))
                                               
    # ASPCF: cut each chromosome where LRR or BAF changes level
    parts <- mclapply(seq_along(samples), function(i) {
      part_file <- file.path(work, paste0(samples[i], ".segmented.rds"))   # kept, so a stopped run resumes
      if (file.exists(part_file)) return(readRDS(part_file))
      one <- ascat.aspcf(loaded, selectsamples = i, ascat.gg = germline, out.dir = work, seed = SEED)
      # keep this array's own column only: ascat.aspcf returns a table as wide as all arrays,
      # and 47 of those sent back to the main process together do not fit in memory
      part <- list(logr = one$Tumor_LogR_segmented[, i], baf = one$Tumor_BAF_segmented[[i]])
      saveRDS(part, part_file)
      part
    }, mc.cores = cores)
    failed <- !vapply(parts, is.list, logical(1))
    if (any(failed)) stop("segmentation failed for ", toString(samples[failed]), ": ", parts[[which(failed)[1]]])
    # put the arrays together again
    segmented <- loaded
    segmented$Tumor_LogR_segmented <- vapply(parts, `[[`, numeric(nrow(loaded$Tumor_LogR)), "logr")
    dimnames(segmented$Tumor_LogR_segmented) <- dimnames(loaded$Tumor_LogR)
    segmented$Tumor_BAF_segmented <- lapply(parts, `[[`, "baf")
    segmented$failedarrays <- germline$failedarrays
    saveRDS(segmented, segmented_file)
    message("segmented in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
  }

  settings <- list(default = list(gamma = DEFAULT_GAMMA, fixed = FALSE),
                   fixed = list(gamma = DEFAULT_GAMMA, fixed = TRUE),
                   calibrated = list(gamma = gamma_x, fixed = TRUE))
  fits <- list()
  for (name in names(settings)) {
    plot_dir <- file.path(work, "plots", name)
    dir.create(plot_dir, showWarnings = FALSE)
    # one after the other (about 1 s per array): ASCAT draws plots here, and drawing
    # in parallel workers (mclapply) hangs
    res <- lapply(seq_along(samples), function(i) {
      run_one(segmented, i, settings[[name]]$gamma, settings[[name]]$fixed, plot_dir)
    })
    segments <- rbindlist(lapply(res, `[[`, "segments"))
    fits[[name]] <- cbind(setting = name, rbindlist(lapply(res, `[[`, "fit")))
    fwrite(segments, file.path(out_dir, sprintf("segments_ASCAT_%s.tsv", name)), sep = "\t")
    message(sprintf("%-10s %d of %d arrays with a solution, %d segments", name,
                    sum(fits[[name]]$status != "no solution"), length(samples), nrow(segments)))
  }
  fwrite(rbindlist(fits), file.path(out_dir, "fit_ASCAT.tsv"), sep = "\t")
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
