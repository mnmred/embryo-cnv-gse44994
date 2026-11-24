# Step 07, evaluation figures - the figures of the talk that need the pipeline's own files.
# The ones here need what only the pipeline writes: CEL
# intensities, probe-level files, ASPCF segments, noise and BAF QC.
#
#   P1_table2_embryo_classes    the paper's Table 2 for the deposited embryos     (slide: dataset)
#   01_raw_intensity            probe intensities before normalisation           (slide: normalise)
#   02_normalised_intensity     probe intensities after quantile normalisation   (slide: normalise)
#   04_lrr_noise_qc             SNP noise (MAPD) against megabase-scale waves    (slide: noise)
#   05_wave_correction          wave SD before and after the GC correction       (slide: waves)
#   08_genotype_clusters        A and B signal of one SNP, by genotype call      (slide: baf)
#   19_calling_thresholds       segment size against copies, with the thresholds (slide: calling)
#
# Run (repo root, bioenv active):
#   Rscript scripts/R/07_evaluation_figures.R [out_dir] [fig_dir]
#   defaults: data/processed  results/figures/evaluation/pipeline
# Reads   data/raw/cel/ (figure 01), data/meta/GSE44994_series_matrix.txt (genotype calls)
#         out_dir: metadata.tsv, arrays.tsv, call_noise.tsv, baf_qc.tsv, baf.csv, probe_level/,
#                  segments_ASPCF.tsv
# File names must match snparray/paths.py.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

CEL_DIR <- "data/raw/cel"
SERIES_MATRIX <- "data/meta/GSE44994_series_matrix.txt"
METHOD <- "ASPCF"                                # segments_<METHOD>.tsv

SLIDE_W <- 40 / 3; SLIDE_H <- 7.5                # inches, PowerPoint 16:9
BASE <- 16                                       # font size in points
COMPLEX_FROM <- 5                                # changed chromosomes that count as complex aneuploidy
MIN_CHANGE <- 0.4; FULL_CHANGE <- 0.8            # the thresholds of calling.py, for the guide lines
MIN_SIZE_FULL <- 2.36; MIN_SIZE_MOSAIC <- 10     # Mb

INK <- "#0b0b0b"; INK_2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; TRACK <- "#f0efec"
POINT_LIGHT <- "#b5b3ab"
GENOTYPE_COLOURS <- c(AA = "#1baf7a", AB = "#eda100", BB = "#4a3aa7")
VIOLET <- "#4a3aa7"; ORANGE <- "#eb6834"
EMBRYO_CLASSES <- c("Normal or balanced", "Translocation-related imbalance", "Translocation-unrelated imbalance")
# Xiong et al. 2014, Table 2. It prints 261 related reciprocal imbalances, but its rows sum to 260,
# which is also what makes the total of 499 and the 52.10% of the GEO summary.
PAPER_TABLE2 <- list(Robertsonian = c(126, 51, 41), Reciprocal = c(177, 260, 62))


# helpers

chr_name <- function(chr) ifelse(chr == 23, "X", as.character(chr))

percent <- function(x, digits = 0) sprintf("%.*f%%", digits, 100 * x)


# inputs

# SNPs of chromosomes 1-22 and X (= 23), sorted as in lrr.tsv.gz
read_snps <- function(out_dir) {
  snps <- fread(file.path(out_dir, "probe_level", "snp_annotation.tsv"),
                select = c("fsetid", "man_fsetid", "chrom", "physical_pos", "cytoband"),
                colClasses = c(chrom = "character"))
  snps <- snps[chrom %in% c(as.character(1:22), "X")]
  snps[, chr := as.integer(ifelse(chrom == "X", "23", chrom))]
  snps <- snps[, .(fsetid, probe = man_fsetid, chr, pos = physical_pos, band = cytoband)]
  setorder(snps, chr, pos)
  snps[]
}

# bp range covered by the SNPs of a cytoband ("q32" includes q32.1, q32.2 ...)
band_range <- function(on_chr, band_name, snps) {
  here <- snps[chr == on_chr]
  if (band_name == "pter") return(rep(min(here$pos), 2))
  if (band_name == "qter") return(rep(max(here$pos), 2))
  if (band_name %in% c("p10", "q10")) return(rep(min(here$pos[startsWith(here$band, "q")]), 2))
  # no SNPs on the short arm (13, 14, 15, 22): such a band lies at the start of the array's data
  if (startsWith(band_name, "p") && !any(startsWith(here$band, "p"))) return(rep(min(here$pos), 2))

  wanted <- band_name
  repeat {
    hit <- if (grepl(".", wanted, fixed = TRUE) || nchar(wanted) <= 2) {
      startsWith(here$band, wanted)
    } else {
      here$band == wanted | startsWith(here$band, paste0(wanted, "."))
    }
    if (any(hit)) break
    if (nchar(wanted) <= 1) stop("cannot place band ", chr_name(on_chr), band_name)
    # band not on the array or not in hg19: use the next coarser band
    wanted <- if (grepl(".", wanted, fixed = TRUE)) sub("\\..*$", "", wanted) else substr(wanted, 1, nchar(wanted) - 1)
  }
  range(here$pos[hit])
}

# "46,XX,+9(pter-q33),-12(pter-q23)" -> one row per changed region
parse_reported <- function(meta, snps) {
  one_embryo <- function(s, text) {
    body <- sub("^([0-9]+,)?(XX|XY)", "", text)
    tokens <- regmatches(body, gregexpr("[+-]([0-9]{1,2}|X)(\\([^)]*\\)|[pq])?", body))[[1]]
    bare <- gsub("\\([^)]*\\)", "", body)
    if (length(tokens) != lengths(regmatches(bare, gregexpr("[+-]", bare)))) {
      stop("cannot parse the result of ", s, ": ", text)
    }
    rbindlist(lapply(unique(tokens), function(token) {
      name <- sub("^[+-]([0-9]{1,2}|X).*$", "\\1", token)
      on_chr <- if (name == "X") 23L else as.integer(name)
      where <- gsub("[()]", "", substring(token, 2 + nchar(name)))
      where <- sub("([0-9])([pq]ter)", "\\1-\\2", where)    # GEO has one "q32.1qter"
      here <- snps[chr == on_chr]
      bp <- if (where == "") {
        range(here$pos)                                     # whole chromosome
      } else if (where %in% c("p", "q")) {
        range(here$pos[startsWith(here$band, where)])       # whole arm
      } else {
        range(unlist(lapply(strsplit(where, "-", fixed = TRUE)[[1]],
                            function(b) band_range(on_chr, b, snps))))
      }
      direction <- if (startsWith(token, "+")) 1L else -1L
      data.table(sample = s, chr = on_chr, start = bp[1], end = bp[2], direction = direction,
                 whole = where == "", copies = 2L + direction * sum(tokens == token))   # "+21,+21" = 4 copies
    }))
  }
  rbindlist(Map(one_embryo, meta$sample, meta$reported_result))
}

# the parent's translocation: one row per embryo and translocation chromosome
parse_parents <- function(meta, snps) {
  pattern <- "(t|der|rob)\\(([0-9]+);([0-9]+)\\)\\(([pq][0-9.]+);([pq][0-9.]+)"
  out <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
    m <- regmatches(meta$donor_karyotype[i], regexec(pattern, meta$donor_karyotype[i]))[[1]]
    if (length(m) == 0) stop("cannot parse the parent's karyotype of ", meta$sample[i])
    chrs <- as.integer(m[3:4]); bands <- m[5:6]
    # from the breakpoints: GEO's donor_class labels some Robertsonian carriers "reciprocal"
    robertsonian <- all(bands == "q10")
    data.table(sample = meta$sample[i], chr = chrs, band = bands, robertsonian = robertsonian,
               parent = sprintf("%s(%d;%d)", if (robertsonian) "rob" else "t", chrs[1], chrs[2]))
  }))
  where <- Map(band_range, out$chr, out$band, MoreArgs = list(snps = snps))
  out[, `:=`(band_start = vapply(where, `[`, numeric(1), 1), band_end = vapply(where, `[`, numeric(1), 2))]
  out[]
}

read_segments <- function(out_dir, method) {
  seg <- fread(file.path(out_dir, sprintf("segments_%s.tsv", method)))
  seg[, reported := as.logical(reported)]            # pandas writes True / False
  seg[, direction := as.integer(sign(change))]
  seg[, state := fcase(!reported, "Normal",
                       kind == "mosaic" & direction > 0, "Mosaic gain",
                       kind == "mosaic", "Mosaic loss",
                       direction < 0, "Loss",
                       round(abs(change)) >= 2, "Gain of 2+ copies",
                       default = "Gain")]
  seg[]
}

# the authors' genotype calls: SNP x embryo, 0 = AA, 1 = AB, 2 = BB, NA = no call
read_genotypes <- function(meta, snps) {
  calls <- fread(SERIES_MATRIX, skip = "ID_REF", fill = TRUE, showProgress = FALSE)
  row <- match(snps$probe, calls$ID_REF)
  out <- vapply(meta$gsm, function(g) match(calls[[g]][row], c("AA", "AB", "BB")) - 1L, integer(nrow(snps)))
  colnames(out) <- meta$sample
  out
}

# a SNP x embryo table of steps 04 / 05, as a matrix in the order of `snps`
read_matrix <- function(file, meta, snps) {
  x <- fread(file, select = c("probe", meta$sample), showProgress = FALSE)
  if (!identical(x$probe, snps$probe)) stop(file, " is not in the SNP order of snp_annotation.tsv")
  as.matrix(x[, meta$sample, with = FALSE])
}

# one SNP with all three genotypes whose clusters are of median tightness (so not a hand-picked one)
typical_snp <- function(d) {
  g <- d$genotype[d$autosomal, !d$mismatched, drop = FALSE]
  b <- d$baf[d$autosomal, !d$mismatched, drop = FALSE]
  ss <- 0; n <- 0; enough <- TRUE
  for (k in 0:2) {
    m <- !is.na(g) & g == k
    count <- rowSums(m)
    v <- b; v[!m] <- 0
    ss <- ss + rowSums(v^2) - rowSums(v)^2 / pmax(count, 1)
    n <- n + count
    enough <- enough & count >= 8
  }
  spread <- ifelse(enough, ss / n, NA)
  which(d$autosomal)[which.min(abs(spread - median(spread, na.rm = TRUE)))]
}

# intensity quantiles of every array before and after normalisation, and one SNP's probes
read_probe_level <- function(out_dir, meta, fsetid_wanted) {
  dir <- file.path(out_dir, "probe_level")
  probes <- fread(file.path(dir, "pm_probes.tsv"), select = c("fid", "fsetid", "allele"))
  cels <- readLines(file.path(dir, "cel_order.txt"))
  samples <- meta$sample[match(sub("_.*$", "", cels), meta$gsm)]
  rows <- which(probes$fsetid == fsetid_wanted)
  at <- c(low = 0.05, q1 = 0.25, middle = 0.5, q3 = 0.75, high = 0.95)
  summarise <- function(v, j) cbind(data.table(sample = samples[j]), as.data.table(as.list(
    setNames(quantile(v, at, names = FALSE), names(at)))))

  connection <- file(file.path(dir, "pm_log2.f32"), "rb")
  on.exit(close(connection))
  raw <- normalised <- vector("list", length(cels))
  picked <- matrix(NA_real_, length(rows), length(cels), dimnames = list(NULL, samples))
  for (j in seq_along(cels)) {
    after <- readBin(connection, "numeric", n = nrow(probes), size = 4)
    normalised[[j]] <- summarise(after, j)
    picked[, j] <- after[rows]
    before <- affxparser::readCelIntensities(file.path(CEL_DIR, cels[j]), indices = probes$fid)
    raw[[j]] <- summarise(log2(as.vector(before)), j)
  }
  # one value per allele and array, by median polish as in step 03
  allele_signal <- function(a) {
    fit <- suppressWarnings(medpolish(picked[probes$allele[rows] == a, , drop = FALSE],
                                      eps = 0.01, maxiter = 10, trace.iter = FALSE))
    fit$overall + fit$col
  }
  list(raw = rbindlist(raw), normalised = rbindlist(normalised),
       cluster = data.table(sample = samples, A = allele_signal(0), B = allele_signal(1)))
}


# against the authors

# best reciprocal overlap of every row of x with a row of y (same embryo, chr, direction)
best_overlap <- function(x, y) {
  vapply(seq_len(nrow(x)), function(i) {
    k <- which(y$sample == x$sample[i] & y$chr == x$chr[i] & y$direction == x$direction[i])
    shared <- pmin(x$end[i], y$end[k]) - pmax(x$start[i], y$start[k]) + 1
    longest <- pmax(x$end[i] - x$start[i], y$end[k] - y$start[k]) + 1
    max(0, shared / longest)
  }, numeric(1))
}

# truth and calls of one method, each marked with whether the other side has it
compare <- function(truth, segments) {
  truth <- copy(truth); calls <- segments[reported == TRUE]
  truth[, overlap := best_overlap(truth, calls)]
  truth[, found := overlap > 0]
  calls[, overlap := best_overlap(calls, truth)]
  calls[, confirmed := overlap > 0]
  list(truth = truth, calls = calls)
}

load_data <- function(out_dir, method) {
  meta <- fread(file.path(out_dir, "metadata.tsv"), colClasses = c(couple = "character"), encoding = "UTF-8")
  by_sample <- function(file) {
    x <- fread(file.path(out_dir, file))
    setnames(x, 1, "sample")                         # pandas leaves some index columns unnamed
    x[match(meta$sample, sample)]
  }
  d <- list(out_dir = out_dir, meta = meta, method = method)
  d$arrays <- by_sample("arrays.tsv")
  d$noise <- by_sample("call_noise.tsv")
  d$baf_qc <- by_sample("baf_qc.tsv")
  d$snps <- read_snps(out_dir)
  d$autosomal <- d$snps$chr <= 22
  d$mismatched <- as.logical(d$baf_qc$calls_do_not_match_cel)
  d$parents <- parse_parents(meta, d$snps)
  d$segments <- read_segments(out_dir, method)
  if (!setequal(d$segments$sample, meta$sample)) stop("segments and metadata list different embryos")
  d[c("truth", "calls")] <- compare(parse_reported(meta, d$snps), d$segments)

  message("reading BAF and genotype calls")
  d$baf <- read_matrix(file.path(out_dir, "baf.csv"), meta, d$snps)         # not clipped to 0-1
  d$genotype <- read_genotypes(meta, d$snps)
  d$typical <- typical_snp(d)
  message("reading probe intensities and CEL files")
  d$probe_level <- read_probe_level(out_dir, meta, d$snps$fsetid[d$typical])
  d
}


# style

theme_slide <- function(grid = "") {
  line <- element_line(colour = GRID, linewidth = 0.4)
  theme_minimal(base_size = BASE) +
    theme(text = element_text(colour = INK),
          axis.text = element_text(colour = INK_2, size = rel(0.9)),
          axis.title = element_text(colour = INK_2, size = rel(0.95)),
          panel.grid = element_blank(),
          panel.grid.major.x = if (grepl("x", grid)) line else element_blank(),
          panel.grid.major.y = if (grepl("y", grid)) line else element_blank(),
          plot.title = element_text(face = "bold", size = rel(1.3)),
          plot.subtitle = element_text(colour = INK_2, size = rel(0.92), margin = margin(b = BASE * 0.8)),
          plot.caption = element_text(colour = INK_2, size = rel(0.72), hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          strip.text = element_text(colour = INK, face = "bold", size = rel(0.9), hjust = 0),
          legend.title = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = rel(0.88)),
          plot.margin = margin(BASE, BASE * 1.2, BASE * 0.7, BASE))
}

embryo_axis <- function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.62)))
}


# as in the paper

# Table 2 of the paper, for the deposited embryos
fig_table2 <- function(d) {
  q <- unique(d$parents[, .(sample, robertsonian)])
  q[, type := ifelse(robertsonian, "Robertsonian", "Reciprocal")]
  class_geo <- d$meta$reported_class[match(q$sample, d$meta$sample)]
  q[, authors := fcase(grepl("^normal", class_geo), EMBRYO_CLASSES[1],
                       grepl("^Unbalanced translocation", class_geo), EMBRYO_CLASSES[2],
                       default = EMBRYO_CLASSES[3])]
  q[, pipeline := vapply(sample, function(s) {
    changed <- unique(d$calls$chr[d$calls$sample == s])
    if (length(changed) == 0) return(EMBRYO_CLASSES[1])
    if (length(changed) >= COMPLEX_FROM) return(EMBRYO_CLASSES[3])            # complex, as in Table 2
    if (any(changed %in% d$parents$chr[d$parents$sample == s])) EMBRYO_CLASSES[2] else EMBRYO_CLASSES[3]
  }, character(1))]
  sources <- c("Paper, all embryos", "Authors, these embryos", "This pipeline")
  cells <- rbindlist(lapply(names(PAPER_TABLE2), function(t_) {
    here <- q[type == t_]
    count <- function(v) as.integer(table(factor(v, levels = EMBRYO_CLASSES)))
    data.table(type = t_, class = rep(EMBRYO_CLASSES, 3), source = rep(sources, each = 3),
               n = c(PAPER_TABLE2[[t_]], count(here$authors), count(here$pipeline)))
  }))
  cells[, share := n / sum(n), by = .(type, source)]
  cells[, total := sum(n), by = .(type, source)]
  cells[, `:=`(type = factor(type, levels = names(PAPER_TABLE2)),
               class = factor(sub(" imbalance", "\nimbalance", class), levels = rev(sub(" imbalance", "\nimbalance", EMBRYO_CLASSES))),
               source = factor(source, levels = sources))]
  p <- ggplot(cells, aes(source, class, fill = share)) +
    geom_tile(colour = "white", linewidth = 3) +
    geom_text(aes(label = sprintf("%d of %d  (%s)", n, total, percent(share)), colour = share > 0.5),
              size = 17, size.unit = "pt") +
    facet_wrap(vars(type), ncol = 1) +
    scale_fill_gradient(low = TRACK, high = "#3a3936", limits = c(0, 1), guide = "none") +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = INK), guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         caption = sprintf(paste("Paper = Xiong et al. 2014, Table 2. Carrier type from the parent's breakpoints.",
                                 "%d or more changed chromosomes = complex, counted as unrelated."), COMPLEX_FROM)) +
    theme_slide()
  list(plot = p, seen = sprintf(paste("The pipeline puts %d of %d embryos in the authors' class. The deposited set is easier",
                                      "than the study: %s are normal, against %s in the paper."),
                                sum(q$authors == q$pipeline), nrow(q), percent(mean(q$authors == EMBRYO_CLASSES[1])),
                                percent(sum(vapply(PAPER_TABLE2, `[`, numeric(1), 1)) / sum(unlist(PAPER_TABLE2)))))
}


# signal, LRR and BAF

box_per_embryo <- function(d, q, y_title, limits = NULL) {
  q <- copy(q)[, sample := factor(sample, levels = d$meta$sample)]
  ggplot(q, aes(sample)) +
    geom_boxplot(aes(ymin = low, lower = q1, middle = middle, upper = q3, ymax = high),
                 stat = "identity", fill = TRACK, colour = INK_2, width = 0.65, linewidth = 0.4) +
    coord_cartesian(ylim = limits) +
    labs(x = NULL, y = y_title, caption = "Box = middle 50% of the values, line = median, whiskers = 5th to 95th percentile.") +
    theme_slide("y") + embryo_axis()
}

intensity_limits <- function(d) {
  q <- rbind(d$probe_level$raw, d$probe_level$normalised)
  c(min(q$low), max(q$high))
}

fig_raw_intensity <- function(d) {
  q <- d$probe_level$raw
  list(plot = box_per_embryo(d, q, "log2 raw probe intensity", intensity_limits(d)),
       seen = sprintf(paste("The median brightness ranges from %.2f to %.2f log2 units, a %.1f-fold difference",
                            "between arrays, so they cannot be compared directly."),
                      min(q$middle), max(q$middle), 2^(max(q$middle) - min(q$middle))))
}

fig_normalised_intensity <- function(d) {
  q <- d$probe_level$normalised
  list(plot = box_per_embryo(d, q, "log2 normalised probe intensity", intensity_limits(d)),
       seen = sprintf("After quantile normalisation the array medians differ by at most %.3f log2 units.",
                      max(q$middle) - min(q$middle)))
}

fig_noise_qc <- function(d) {
  q <- d$arrays[, .(sample, mapd, wave)]
  q[, class := ifelse(sample %in% d$truth$sample, "Abnormal in GEO", "Normal in GEO")]
  noisiest <- unique(c(q[order(-mapd), head(sample, 3)], q[order(-wave), head(sample, 2)],
                       d$noise[which.max(tau_arm), sample]))
  p <- ggplot(q, aes(mapd, wave)) +
    geom_point(aes(shape = class), size = 3.5, colour = INK_2, stroke = 1) +
    geom_text(data = q[sample %in% noisiest], aes(label = sample), hjust = -0.12, colour = INK,
              size = 12, size.unit = "pt", check_overlap = TRUE) +
    scale_shape_manual(values = c("Abnormal in GEO" = 16, "Normal in GEO" = 1)) +
    scale_x_continuous(expand = expansion(c(0.05, 0.16))) +
    labs(x = "SNP-to-SNP noise (MAPD)", y = "Megabase-scale waves (SD of LRR)",
         caption = "One point per embryo. A wave SD of 0 means below what can be measured. The noisiest embryos are labelled.") +
    theme_slide("xy")
  list(plot = p, seen = sprintf("MAPD is %.2f to %.2f (median %.2f). The array with the strongest SNP noise is %s.",
                                min(q$mapd), max(q$mapd), median(q$mapd), q[which.max(mapd), sample]))
}

fig_wave_correction <- function(d) {
  q <- d$arrays[, .(sample, before = wave_before, after = wave)][order(-before)]
  q[, sample := factor(sample, levels = sample)]
  long <- melt(q, id.vars = "sample", variable.name = "stage")
  long[, stage := factor(ifelse(stage == "before", "Before correction", "After correction"),
                         levels = c("Before correction", "After correction"))]
  p <- ggplot(q, aes(sample)) +
    geom_segment(aes(xend = sample, y = before, yend = after), colour = POINT_LIGHT, linewidth = 1) +
    geom_point(data = long, aes(y = value, shape = stage), colour = INK, size = 2.8, stroke = 0.9) +
    scale_shape_manual(values = c("Before correction" = 1, "After correction" = 16)) +
    labs(x = NULL, y = "Megabase-scale waves (SD of LRR)",
         caption = "Correction for regional GC, PCR fragment length and probe GC (step 04). 0 = below what can be measured.") +
    theme_slide("y") + embryo_axis()
  list(plot = p, seen = sprintf("The median wave SD falls from %.3f to %.3f, and is lower afterwards in %d of %d embryos.",
                                median(q$before), median(q$after), sum(q$after < q$before), nrow(q)))
}

fig_clusters <- function(d) {
  q <- copy(d$probe_level$cluster)
  i <- match(q$sample, d$meta$sample)
  call <- c("AA", "AB", "BB")[d$genotype[d$typical, i] + 1]
  call[is.na(call) | d$mismatched[i]] <- "No usable call"
  q[, genotype := call]
  p <- ggplot(q, aes(A, B)) +
    geom_abline(slope = 1, colour = GRID) +
    geom_point(aes(colour = genotype, shape = genotype), size = 4, stroke = 1.1) +
    scale_colour_manual(values = c(GENOTYPE_COLOURS, "No usable call" = MUTED)) +
    scale_shape_manual(values = c(AA = 16, AB = 17, BB = 15, "No usable call" = 4)) +
    labs(x = "log2 signal of allele A", y = "log2 signal of allele B",
         caption = "One point per embryo. Colour = the authors' genotype call. SNP chosen automatically as one of median cluster tightness.") +
    theme_slide("xy") + theme(legend.position = "right")
  n <- table(factor(call, levels = c("AA", "AB", "BB")))
  list(plot = p, seen = sprintf("A typical SNP, %s: %d AA, %d AB and %d BB embryos. The clusters are where BAF 0, 0.5 and 1 are anchored.",
                                d$snps$probe[d$typical], n[["AA"]], n[["AB"]], n[["BB"]]))
}


# segments to calls

fig_thresholds <- function(d) {
  q <- d$segments[size_mb >= 1, .(sample, chr, start, end, size_mb, change, reported)]
  q[d$calls, on = .(sample, chr, start, end), confirmed := i.confirmed]
  q[, status := fcase(!reported, "Not called", confirmed, "Called, agrees with the authors",
                      default = "Called, not reported by the authors")]
  setorder(q, reported)
  p <- ggplot(q, aes(size_mb, change)) +
    geom_hline(yintercept = c(-1, 1) * MIN_CHANGE, colour = MUTED, linetype = "22") +
    geom_hline(yintercept = c(-1, 1) * FULL_CHANGE, colour = MUTED, linetype = "13") +
    geom_vline(xintercept = c(MIN_SIZE_FULL, MIN_SIZE_MOSAIC), colour = MUTED, linetype = "42") +
    geom_point(aes(colour = status, shape = status, size = status)) +
    annotate("text", x = min(q$size_mb), y = c(MIN_CHANGE, FULL_CHANGE), label = c("mosaic", "full"),
             hjust = 0, vjust = -0.5, colour = INK_2, size = 13, size.unit = "pt") +
    annotate("text", x = c(MIN_SIZE_FULL, MIN_SIZE_MOSAIC), y = min(q$change), hjust = -0.1,
             label = sprintf("%g Mb", c(MIN_SIZE_FULL, MIN_SIZE_MOSAIC)), colour = INK_2, size = 13, size.unit = "pt") +
    scale_colour_manual(values = c("Not called" = POINT_LIGHT, "Called, agrees with the authors" = VIOLET,
                                   "Called, not reported by the authors" = ORANGE)) +
    scale_shape_manual(values = c("Not called" = 16, "Called, agrees with the authors" = 16,
                                  "Called, not reported by the authors" = 17)) +
    scale_size_manual(values = c("Not called" = 1.2, "Called, agrees with the authors" = 3.3,
                                 "Called, not reported by the authors" = 3.6)) +
    scale_x_log10(breaks = c(1, 3, 10, 30, 100, 250), labels = c("1", "3", "10", "30", "100", "250")) +
    labs(x = "Segment size (Mb, log scale)", y = "Copies gained (+) or lost (-)",
         caption = "One point per segment of at least 1 Mb. A call also needs |z| ≥ 3, so some points beyond the lines are not called.") +
    theme_slide()
  extra <- d$calls[confirmed == FALSE]
  weak <- MIN_CHANGE + 0.15                        # "close to the threshold"
  if (nrow(extra) == 0) {
    below <- sum(!d$truth$found)
    return(list(plot = p, seen = sprintf(
      "All %d calls are also reported by the authors. The price is recall: %d reported changes are not called.",
      nrow(d$calls), below)))
  }
  list(plot = p, seen = sprintf("The authors do not report %d of the calls: %d gains and %d losses. %d of them change by less than %.2f copies, close to the %.1f threshold.",
                                nrow(extra), sum(extra$direction > 0), sum(extra$direction < 0),
                                sum(abs(extra$change) < weak), weak, MIN_CHANGE))
}


# output

save_slide <- function(plot, target, dpi = 200) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(target, SLIDE_W, SLIDE_H, units = "in", res = dpi, background = "white")
  on.exit(dev.off())
  print(plot)
}

wrap <- function(text, width = 118) paste(strwrap(text, width = width), collapse = "\n")


# what to draw

FIGURES <- list(
  list(id = "P1", name = "table2_embryo_classes", title = "Paper's Table 2: classes of embryos", make = fig_table2),
  list(id = "01", name = "raw_intensity", title = "Raw probe intensities", make = fig_raw_intensity),
  list(id = "02", name = "normalised_intensity", title = "Probe intensities after quantile normalisation",
       make = fig_normalised_intensity),
  list(id = "04", name = "lrr_noise_qc", title = "Quality control: LRR noise per embryo", make = fig_noise_qc),
  list(id = "05", name = "wave_correction", title = "Amplification waves before and after correction",
       make = fig_wave_correction),
  list(id = "08", name = "genotype_clusters", title = "Genotype clusters of one SNP", make = fig_clusters),
  list(id = "19", name = "calling_thresholds", title = "From segments to calls: size and effect thresholds",
       make = fig_thresholds)
)

main <- function(args) {
  out_dir <- if (length(args) >= 1) args[1] else "data/processed"
  fig_dir <- if (length(args) >= 2) args[2] else "results/figures/evaluation/pipeline"
  d <- load_data(out_dir, METHOD)
  for (f in FIGURES) {
    out <- f$make(d)
    target <- file.path(fig_dir, sprintf("%s_%s.png", f$id, f$name))
    save_slide(out$plot + labs(title = f$title, subtitle = wrap(out$seen)), target)
    message(sprintf("%-32s %s", basename(target), out$seen))
  }
  message(length(FIGURES), " figures -> ", fig_dir)
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
