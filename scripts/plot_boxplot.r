#!/usr/bin/env Rscript
# Count targets/decoys passing FDR (q-value) thresholds in results/3_pep outputs.
# Base R only.

# -------- args --------
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k, d=NULL){ i <- which(args==k); if(!length(i)||i==length(args)) d else args[i+1] }

indir  <- arg("--indir",  "results/3_pep")
outdir <- arg("--outdir", file.path(indir, "results/summary_boxplot"))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
seed_in <- arg("--seed", NA)

# Thresholds to report
thr <- c(0.01, 0.05)

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# -------- helpers --------
seed_from_name <- function(x){
  # matches ..._seed1_... or ...seed1...
  b <- basename(x)
  m <- regexpr("seed([0-9]+)", b, perl=TRUE)
  if (m < 0) NA_integer_ else as.integer(sub(".*seed([0-9]+).*","\\1", b))
}

read_tsv_loose <- function(p){
  con <- file(p, open="r"); on.exit(close(con))
  header <- strsplit(readLines(con, n=1, warn=FALSE), "\t", fixed=TRUE)[[1]]
  dat <- tryCatch(
    read.table(p, sep="\t", header=TRUE, quote="", comment.char="",
               check.names=FALSE, stringsAsFactors=FALSE),
    error=function(e) stop("Failed to read ", basename(p), ": ", e$message)
  )
  if (length(names(dat)) != length(header)) names(dat) <- header
  dat
}

norm_name <- function(x){
  tolower(gsub("[[:space:]_\\-\\.\\(\\)]+", "", x))
}

find_fdr_col <- function(cols){
  n <- norm_name(cols)
  # preferred exact matches
  wanted <- c("qvalue", "qval", "q", "fdr", "qvalues", "qvalueestimated")
  for (w in wanted){
    hit <- which(n == w)
    if (length(hit)) return(cols[hit[1]])
  }
  # fallback: anything containing "qvalue" or ending with "fdr" / "q"
  hit <- grep("qvalue", n)
  if (length(hit)) return(cols[hit[1]])
  hit <- grep("fdr$", n)
  if (length(hit)) return(cols[hit[1]])
  NA_character_
}

find_label_col <- function(cols){
  n <- norm_name(cols)
  hit <- which(n == "label")
  if (length(hit)) return(cols[hit[1]])
  NA_character_
}

# -------- discover files under results/3_pep --------
# Accept both spellings: ipspline / ispline (your earlier scripts check both) :contentReference[oaicite:1]{index=1}
files <- unique(c(
  Sys.glob(file.path(indir, "*_pep_seed*_ipspline.txt")),
  Sys.glob(file.path(indir, "*_pep_seed*_ispline.txt")),
  Sys.glob(file.path(indir, "*_pep_ipspline.txt")),
  Sys.glob(file.path(indir, "*_pep_ispline.txt"))
))

if (!length(files)) {
  stop("No pyIsoPEP outputs found under ", indir, " (expected *_pep*_ipspline.txt or *_pep*_ispline.txt).", call.=FALSE)
}

# optional seed filter
if (!is.na(seed_in)) {
  s <- as.integer(seed_in)
  files <- files[seed_from_name(files) == s]
  if (!length(files)) stop("No files found for seed ", s, " under ", indir, call.=FALSE)
}

# -------- count per file (targets + decoys) --------
rows <- list()

for (f in files){
  df <- read_tsv_loose(f)

  seed <- seed_from_name(f)
  if (is.na(seed)) seed <- -1L  # if filenames don’t encode seed

  fdr_col <- find_fdr_col(names(df))
  if (is.na(fdr_col)) {
    stop("No FDR/q-value column found in ", basename(f),
         ". Looked for variants of: q, qval, qvalue, fdr.", call.=FALSE)
  }

  label_col <- find_label_col(names(df))
  if (is.na(label_col)) {
    stop("No Label column found in ", basename(f), ". Needed to split target vs decoy.", call.=FALSE)
  }

  q <- suppressWarnings(as.numeric(df[[fdr_col]]))
  lab <- suppressWarnings(as.integer(df[[label_col]]))
  lab[is.na(lab)] <- 0L

  # define groups:
  # target = 1
  # decoy  = -1 or 0 (depending on how upstream wrote it; we treat both as decoy-like)
  is_target <- lab == 1L
  is_decoy  <- lab != 1L

  for (t in thr){
    n_target_total <- sum(!is.na(q) & is_target)
    n_target_pass  <- sum(!is.na(q) & is_target & q <= t)

    n_decoy_total <- sum(!is.na(q) & is_decoy)
    n_decoy_pass  <- sum(!is.na(q) & is_decoy & q <= t)

    rows[[length(rows)+1]] <- data.frame(
      seed=seed,
      file=basename(f),
      fdr_col=fdr_col,
      threshold=t,
      target_pass=n_target_pass,
      target_total=n_target_total,
      decoy_pass=n_decoy_pass,
      decoy_total=n_decoy_total,
      stringsAsFactors=FALSE
    )
  }
}

summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$seed, summary_df$threshold, summary_df$file), ]

# Write CSV summary
csv_path <- file.path(outdir, "fdr_pass_summary.csv")
write.table(summary_df, file=csv_path, sep=",", row.names=FALSE)
cat("✅ wrote:", csv_path, "\n")

# Aggregate per seed across files
agg <- aggregate(cbind(target_pass, target_total, decoy_pass, decoy_total) ~ seed + threshold,
                 data=summary_df, FUN=sum, na.rm=TRUE)
agg_path <- file.path(outdir, "fdr_pass_summary_by_seed.csv")
write.table(agg, file=agg_path, sep=",", row.names=FALSE)
cat("✅ wrote:", agg_path, "\n")

# ---- ggplot2 boxplot across seeds, with per-seed points ----
suppressPackageStartupMessages({
  library(ggplot2)
})

# Convert aggregated per-seed counts into long format (base R)
# agg has columns: seed, threshold, target_pass, target_total, decoy_pass, decoy_total
long <-  data.frame(seed=agg$seed, threshold=agg$threshold, group="targets", pass=agg$target_pass, stringsAsFactors=FALSE)

long$threshold <- factor(
  long$threshold,
  levels=c(0.01, 0.05),
  labels=c("FDR ≤ 1%", "FDR ≤ 5%")
)
# long$group <- factor(long$group, levels=c("targets","decoys"))

p <- ggplot(long, aes(x = threshold, y = pass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    aes(group = seed),
    alpha = 0.5,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(after_stat(y), 0)),
    vjust = -1,
    size = 4
  ) +
  theme_bw() +
  labs(
    title = "Target counts passing FDR thresholds (across seeds)",
    x = NULL,
    y = "Target count passing threshold"
  )

out_png <- file.path(outdir, "fdr_pass_boxplot.png")
ggsave(out_png, p, width=8.5, height=5.0, dpi=150)
cat("✅ wrote:", out_png, "\n")
