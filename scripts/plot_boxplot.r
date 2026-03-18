#!/usr/bin/env Rscript
# Count peptides passing PEP and q-value thresholds in pyIsoPEP outputs.
# Reports counts per seed and per file-basename (everything before 'seed' in filename).

# -------- args --------
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k, d=NULL){
  i <- which(args == k)
  if(!length(i) || i == length(args)) d else args[i+1]
}

# ---- determine project root to avoid creating folders under scripts/ ----
get_script_dir <- function(){
  ca <- commandArgs(trailingOnly=FALSE)
  f <- grep("^--file=", ca, value=TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=","", f[1]))))
  # fallback: current working directory (less ideal)
  getwd()
}

script_dir <- get_script_dir()
project_root <- if (basename(script_dir) == "scripts") dirname(script_dir) else getwd()

indir  <- arg("--indir",  file.path(project_root, "results", "3_pep"))
outdir <- arg("--outdir", file.path(project_root, "results", "summary_boxplot"))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

seed_in <- arg("--seed", NA)

# Thresholds to report
thr <- c(0.01, 0.05)

# Column names requested (exact)
pep_col_name  <- "pyIsoPEP PEP"
qval_col_name <- "pyIsoPEP q-value from FDR"

# -------- helpers --------
seed_from_name <- function(x){
  b <- basename(x)
  m <- regexpr("seed([0-9]+)", b, perl=TRUE)
  if (m < 0) NA_integer_ else as.integer(sub(".*seed([0-9]+).*","\\1", b))
}

# basename = everything before 'seed' in the filename (trim trailing separators)
base_from_name <- function(x){
  b <- basename(x)
  pre <- sub("(.*)seed[0-9]+.*", "\\1", b, perl=TRUE)
  pre <- sub("[_\\-.]+$", "", pre, perl=TRUE)
  if (identical(pre, b)) pre <- sub("\\.[^.]+$", "", b)  # if no seed found, use stem
  pre
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

find_col_exact_or_norm <- function(cols, target_exact){
  # prefer exact match first
  hit <- which(cols == target_exact)
  if (length(hit)) return(cols[hit[1]])

  # fallback: normalized match (helps if upstream altered whitespace a bit)
  n <- norm_name(cols)
  t <- norm_name(target_exact)
  hit <- which(n == t)
  if (length(hit)) return(cols[hit[1]])

  NA_character_
}

find_label_col <- function(cols){
  n <- norm_name(cols)
  hit <- which(n == "label")
  if (length(hit)) return(cols[hit[1]])
  NA_character_
}

# -------- discover files under indir --------
files <- unique(c(
  Sys.glob(file.path(indir, "*_pep_seed*_ipspline.txt")),
  Sys.glob(file.path(indir, "*_pep_seed*_ispline.txt")),
  Sys.glob(file.path(indir, "*_pep_ipspline.txt")),
  Sys.glob(file.path(indir, "*_pep_ispline.txt"))
))

if (!length(files)) {
  stop("No pyIsoPEP outputs found under ", indir,
       " (expected *_pep*_ipspline.txt or *_pep*_ispline.txt).", call.=FALSE)
}

# optional seed filter
if (!is.na(seed_in)) {
  s <- as.integer(seed_in)
  files <- files[seed_from_name(files) == s]
  if (!length(files)) stop("No files found for seed ", s, " under ", indir, call.=FALSE)
}

# -------- count per file-basename + seed --------
rows <- list()

for (f in files){
  df <- read_tsv_loose(f)

  seed <- seed_from_name(f)
  if (is.na(seed)) seed <- -1L

  base <- base_from_name(f)

  pep_col  <- find_col_exact_or_norm(names(df), pep_col_name)
  qval_col <- find_col_exact_or_norm(names(df), qval_col_name)

  if (is.na(pep_col)) {
    stop("Missing column '", pep_col_name, "' in ", basename(f), call.=FALSE)
  }
  if (is.na(qval_col)) {
    stop("Missing column '", qval_col_name, "' in ", basename(f), call.=FALSE)
  }

  # If Label exists, count only targets (Label==1). Otherwise count all rows.
  label_col <- find_label_col(names(df))
  if (!is.na(label_col)) {
    lab <- suppressWarnings(as.integer(df[[label_col]]))
    is_target <- !is.na(lab) & (lab == 1L)
  } else {
    is_target <- rep(TRUE, nrow(df))
  }

  pep  <- suppressWarnings(as.numeric(df[[pep_col]]))
  qval <- suppressWarnings(as.numeric(df[[qval_col]]))

  for (t in thr){
    rows[[length(rows)+1]] <- data.frame(
      seed = seed,
      basename = base,
      file = basename(f),
      threshold = t,

      pep_total = sum(is_target & !is.na(pep)),
      pep_pass  = sum(is_target & !is.na(pep)  & pep  <= t),

      q_total   = sum(is_target & !is.na(qval)),
      q_pass    = sum(is_target & !is.na(qval) & qval <= t),

      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$basename, summary_df$seed, summary_df$threshold, summary_df$file), ]

# Write detailed CSV summary (per file)
csv_path <- file.path(outdir, "pep_q_pass_summary_by_file.csv")
write.table(summary_df, file=csv_path, sep=",", row.names=FALSE)
cat("✅ wrote:", csv_path, "\n")

# Aggregate per (basename, seed) across files
agg <- aggregate(cbind(pep_pass, pep_total, q_pass, q_total) ~ basename + seed + threshold,
                 data=summary_df, FUN=sum, na.rm=TRUE)
agg_path <- file.path(outdir, "pep_q_pass_summary_by_basename_seed.csv")
write.table(agg, file=agg_path, sep=",", row.names=FALSE)
cat("✅ wrote:", agg_path, "\n")

# Optional: quick boxplots across seeds (targets passing thresholds)
# Produces one plot for PEP and one for q-value
suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly=TRUE)
})



if (ok) {
  library(ggplot2)

make_plot <- function(metric_label, pass_col){
  long <- data.frame(
    basename = agg$basename,
    seed = agg$seed,
    threshold = factor(agg$threshold,
                       levels = c(0.01, 0.05),
                       labels = c("≤ 1%", "≤ 5%")),
    pass = agg[[pass_col]],
    stringsAsFactors = FALSE
  )

  facet_order <- c(
    "drep_tcell100p_scored_target",
    "tcell30p_scored_target",
    "tcell10p_scored_target",
    "tcell3p_scored_target",
    "tcell1p_scored_target",
    "tcell0pt3p_scored_target"
  )

  long$basename <- factor(long$basename, levels = facet_order)

  ggplot(long, aes(x = threshold, y = pass)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
      aes(group = seed),
      alpha = 0.5,
      position = position_jitter(width = 0.12, height = 0)
    ) +
    facet_wrap(~ basename, scales = "free_y") + ylim(0,132) + 
    theme_bw() + hline(yintercept = 132, linetype = "dashed", color = "gray") +
    labs(
      title = paste0(metric_label, ": target counts passing threshold (across seeds)"),
      x = NULL,
      y = "Count passing threshold"
    )
}


  p1 <- make_plot("PEP", "pep_pass")
  p1_path <- file.path(outdir, "pep_pass_boxplot_by_basename.png")
  ggsave(p1_path, p1, width=11, height=6, dpi=150)
  cat("✅ wrote:", p1_path, "\n")

  p2 <- make_plot("q-value from FDR", "q_pass")
  p2_path <- file.path(outdir, "q_pass_boxplot_by_basename.png")
  ggsave(p2_path, p2, width=11, height=6, dpi=150)
  cat("✅ wrote:", p2_path, "\n")
} else {
  cat("ℹ️ ggplot2 not installed; skipped PNG boxplots. CSV summaries were still written.\n")
}
