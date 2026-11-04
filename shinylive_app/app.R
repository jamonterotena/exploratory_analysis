# =============================================================================
# Tomato crosses: genomic features vs. recombination -- interactive explorer
# Shinylive version: runs entirely in the visitor's browser (R via webR), so it
# can be hosted as a static website, e.g. on GitHub Pages.
#
# Build the website (in normal R, from the folder that contains this app folder):
#   install.packages("shinylive")
#   shinylive::export("shinylive_app", "docs")
#   httpuv::runStaticServer("docs")   # test locally, then push docs/ to GitHub
#
# Data are never part of the website: each user loads their own BED files, which
# stay in their browser and are not sent anywhere. Do NOT copy .bed files into
# this folder, or they will be published with the site.
# No password: the site's code is public, so a password could not protect it.
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(DT)
})

options(shiny.maxRequestSize = 4 * 1024^3)   # uploads stay in the browser
setDTthreads(1)                              # no multithreading in the browser

# ---- Configuration ----------------------------------------------------------
# Columns that are NOT plottable features. Everything else that is numeric is
# treated as a feature automatically, so differing feature columns between
# files need no code changes.
KEY_COLS    <- c("Chr_ID", "Interval_start", "Interval_end", "Dataset",
                 "Background_parent", "Donor_parent", "Male")   # used to merge files
META_COLS   <- c(KEY_COLS, "Bin_size_bp", "Number_of_COs", "Number_of_hybrids",
                 "Expected_COs", "Obs_exp_ratio", "Poisson_p",
                 "Interval_midpoint", "Spot_class")
TARGET_COLS <- c("cM_Mbp", "cM_Mbp_sex_avg")

# Suffix of the recombination-rate copies restricted to recombinogenic bins
# (bins with >= 1 CO), used in the distributions tab
REC_SUFFIX <- "_recombinogenic"
is_rate_col <- function(v) v %in% c(TARGET_COLS, paste0(TARGET_COLS, REC_SUFFIX))

# Variables plotted by default, in this order (missing ones are skipped).
DEFAULT_PLOT_VARS <- c("cM_Mbp", "SYN", "INV",
                       "GA_P1", "TE_P1", "CHR_P1", "POS_P1",
                       "GA_P2", "TE_P2", "CHR_P2", "POS_P2")

SPOT_COLORS <- c(Hotspot = "red", Coldspot = "blue")

# ggplotly draws the hot/cold spot shading as one trace per facet panel, each
# with its own legend entry. Put all traces of one region type in a single
# legend group (clicking it toggles that type in every panel) and show the
# entry only once, under a clean name.
merge_spot_legend <- function(pl) {
  seen <- character(0)
  for (i in seq_along(pl$x$data)) {
    tr <- pl$x$data[[i]]
    nm <- tr$name %||% ""
    type <- names(SPOT_COLORS)[vapply(names(SPOT_COLORS),
                                      function(s) grepl(s, nm, fixed = TRUE), logical(1))]
    if (length(type) != 1) next
    tr$name        <- type
    tr$legendgroup <- type
    tr$showlegend  <- !(type %in% seen)
    seen <- c(seen, type)
    pl$x$data[[i]] <- tr
  }
  pl
}

# Diverging scale for correlations: blue (negative) - neutral grey (0) - red
# (positive), matching coldspot = blue / hotspot = red.
CORR_POLES       <- c(neg = "#1c5cab", mid = "#f0efec", pos = "#b8302f")
CORR_TABLE_POLES <- c(neg = "#86b6ef", mid = "#f0efec", pos = "#ee9a99")  # lighter, for cell backgrounds

# ---- Helpers ----------------------------------------------------------------

# Read a BED-like file whose leading lines start with '#', followed by a header
read_bed <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  n <- 0L
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l)) stop("No header line found in ", path)
    if (!startsWith(l, "#")) break
    n <- n + 1L
  }
  fread(path, skip = n, header = TRUE, na.strings = c("NA", ""), data.table = TRUE)
}

# Merge several tables on the key columns; columns already present are not duplicated
combine_tables <- function(dts) {
  Reduce(function(a, b) {
    keys <- intersect(KEY_COLS, intersect(names(a), names(b)))
    b <- b[, c(keys, setdiff(names(b), names(a))), with = FALSE]
    if (ncol(b) == length(keys)) return(a)
    merge(a, b, by = keys, all = TRUE, sort = FALSE)
  }, dts)
}

# Assign a feature column to a category from the tokens in its name
# (split on non-alphanumerics, e.g. "TE_Gypsy_P1" -> TE, GYPSY, P1).
# Adjust the patterns below if your column names use different prefixes.
CATEGORY_LEVELS <- c("Structural variants", "Genome", "Transposable elements", "Chromatin", "Repeats",
                     "Centromeres", "Other")
feature_category <- function(x) {
  tk <- strsplit(toupper(x), "[^A-Z0-9]+")[[1]]
  has <- function(re) any(grepl(re, tk))
  if (has("^CEN[0-9]*$"))                          "Centromeres"
  else if (has("^CHR[0-9]*$|COMPART|EIGEN"))       "Chromatin"
  else if (has("^REP$|REPEAT"))                    "Repeats"
  else if (has("^TE[0-9]*$"))                      "Transposable elements"
  else if (has("^GA[0-9]*$"))                      "Genome"
  else if (has("^(SYN|SYNAL|HDR|TDM|INV|NLO)$"))   "Structural variants"
  else                                             "Other"
}

rescale01 <- function(x) {
  if (all(is.na(x))) return(x)
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - r[1]) / diff(r)
}

# Hotspots / coldspots from a single ordered series of sliding-window bins.
# A region is called when the mean Poisson_p over a run of consecutive bins
# spanning >= min_len bp is < hot_p (hotspot) or >= cold_p (coldspot).
# Overlapping qualifying runs are merged into one region.
call_spots <- function(d, min_len, hot_p, cold_p) {
  n <- nrow(d)
  empty <- data.table(start = numeric(), end = numeric(), type = character())
  cls <- rep(NA_character_, n)
  if (n < 2 || !"Poisson_p" %in% names(d)) return(list(regions = empty, class = cls))
  
  step <- median(diff(d$Interval_start), na.rm = TRUE)
  size <- median(d$Interval_end - d$Interval_start, na.rm = TRUE)
  if (!isTRUE(step > 0) || !is.finite(size) || !isTRUE(is.finite(min_len)))
    return(list(regions = empty, class = cls))
  k <- max(1L, as.integer(ceiling((min_len - size) / step)) + 1L)
  if (n < k) return(list(regions = empty, class = cls))
  
  m  <- frollmean(d$Poisson_p, k, align = "left")
  s2 <- shift(d$Interval_start, k - 1L, type = "lead")
  contiguous <- !is.na(s2) & (s2 - d$Interval_start) <= (k - 1L) * step * 1.001
  
  covered <- function(win) {
    f <- rep(FALSE, n)
    for (i in which(win)) f[i:(i + k - 1L)] <- TRUE
    f
  }
  hot  <- covered(contiguous & !is.na(m) & m < hot_p)
  cold <- covered(contiguous & !is.na(m) & m >= cold_p - 1e-9) & !hot
  
  runs <- function(flag, type) {
    r <- rle(flag)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    keep <- r$values
    if (!any(keep)) return(empty)
    data.table(start = d$Interval_start[starts[keep]],
               end   = d$Interval_end[ends[keep]], type = type)
  }
  cls[hot]  <- "Hotspot"
  cls[cold] <- "Coldspot"
  list(regions = rbind(runs(hot, "Hotspot"), runs(cold, "Coldspot")), class = cls)
}

# Rows for one donor / sex / dataset (all chromosomes)
subset_rows <- function(d, sel_donor, sel_male, sel_dataset) {
  d <- d[Donor_parent == sel_donor & Dataset == sel_dataset]
  if (sel_male == "NA") d[is.na(Male)] else d[!is.na(Male) & Male == as.numeric(sel_male)]
}

# Non-overlapping tiling of sliding-window bins (sorted by start): greedily keep
# each bin that starts at or after the end of the last kept bin.
tile_keep <- function(start, end) {
  keep <- logical(length(start)); last_end <- -Inf
  for (i in seq_along(start)) if (isTRUE(start[i] >= last_end) && !is.na(end[i])) {
    keep[i] <- TRUE
    last_end <- end[i]
  }
  keep
}

has_token <- function(x, token) token %in% strsplit(toupper(x), "[^A-Z0-9]+")[[1]]

# Default presence threshold: 0.50 for inversions and chromatin state, else the median
default_threshold <- function(f, x) {
  if (has_token(f, "INV") || feature_category(f) == "Chromatin") 0.5 else median(x, na.rm = TRUE)
}

# ---- Feature vs. CO state test ----------------------------------------------
# Inputs are non-overlapping bins: CO counts, bin lengths (bp), feature presence,
# CO state (TRUE = the state of interest) and chromosome.
# Poisson test: observed COs in bins with the feature vs. the number expected if
# each chromosome's COs were spread evenly along its length (so pooled results
# are stratified by chromosome). 2x2 table of feature presence vs. CO state:
# precision = P(state | feature), recall = P(feature | state), odds ratio
# (0.5 added to every cell when one is zero) with Fisher's exact test.
feature_state_test <- function(co, len, pos, state, chr, rate = NULL) {
  co <- as.numeric(co); len <- as.numeric(len)    # avoid integer overflow in products
  if (is.null(rate)) rate <- rep(NA_real_, length(co))
  ok <- is.finite(co) & is.finite(len) & !is.na(pos) & !is.na(chr)
  co <- co[ok]; len <- len[ok]; pos <- pos[ok]; chr <- chr[ok]; rate <- as.numeric(rate[ok])
  state <- state[ok]; state[is.na(state)] <- FALSE
  Lc <- tapply(len, chr, sum)
  Pc <- tapply(len * pos, chr, sum)
  Oc <- tapply(co, chr, sum)
  E <- sum(Oc * Pc / Lc)
  O <- sum(co[pos])
  n11 <- sum(pos & state);  n10 <- sum(pos & !state)
  n01 <- sum(!pos & state); n00 <- sum(!pos & !state)
  full <- isTRUE((n11 + n10) > 0 && (n01 + n00) > 0 && (n11 + n01) > 0 && (n10 + n00) > 0)
  or <- if (min(n11, n10, n01, n00) == 0)
    (n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)) else n11 * n00 / (n10 * n01)
  hasE <- isTRUE(E > 0)
  data.frame(
    Bins_pos  = n11 + n10, Bins = length(pos), O = O, E = E,
    OE        = if (hasE) O / E else NA_real_,
    Poisson_p = if (hasE) poisson.test(round(O), E)$p.value else NA_real_,
    Precision = if (n11 + n10 > 0) n11 / (n11 + n10) else NA_real_,
    Recall    = if (n11 + n01 > 0) n11 / (n11 + n01) else NA_real_,
    OR        = if (full) or else NA_real_,
    Fisher_p  = if (full) fisher.test(matrix(c(n11, n01, n10, n00), 2))$p.value else NA_real_,
    rate_tests(co, len, pos, chr, rate))
}

# Recombination rate with vs. without the feature.
#  - Mean rate (all bins) in each subset.
#  - Rate ratio and p-value from a quasi-Poisson GLM, COs ~ feature (+ chromosome)
#    with log(bin length) as offset. The quasi-Poisson dispersion absorbs the
#    extra variance from the many zero bins and the hotspots, so this p-value is
#    more conservative than the plain Poisson test.
#  - Mean rate in recombinogenic bins (>= 1 CO) in each subset, with a Wilcoxon
#    rank-sum test: does the feature change how much recombination happens
#    where it happens at all?
rate_tests <- function(co, len, pos, chr, rate) {
  na <- data.frame(Rate_pos = NA_real_, Rate_neg = NA_real_, QP_RR = NA_real_,
                   QP_p = NA_real_, RecRate_pos = NA_real_, RecRate_neg = NA_real_,
                   Wilcox_p = NA_real_)
  if (!any(pos) || all(pos)) return(na)
  out <- na
  out$Rate_pos <- mean(rate[pos], na.rm = TRUE)
  out$Rate_neg <- mean(rate[!pos], na.rm = TRUE)
  use <- len > 0
  if (sum(co[use]) > 0) {
    f <- if (length(unique(chr[use])) > 1) co ~ pos + factor(chr) else co ~ pos
    fit <- tryCatch(suppressWarnings(
      glm(f, family = quasipoisson(), offset = log(len),
          data = data.frame(co = co[use], pos = pos[use], chr = chr[use], len = len[use]))),
      error = function(e) NULL)
    if (!is.null(fit) && "posTRUE" %in% rownames(coef(summary(fit)))) {
      cs <- coef(summary(fit))["posTRUE", ]
      out$QP_RR <- exp(cs[["Estimate"]]); out$QP_p <- cs[[4]]
    }
  }
  rec <- co > 0 & is.finite(rate)
  a <- rate[rec & pos]; b <- rate[rec & !pos]
  out$RecRate_pos <- if (length(a)) mean(a) else NA_real_
  out$RecRate_neg <- if (length(b)) mean(b) else NA_real_
  if (length(a) >= 2 && length(b) >= 2)
    out$Wilcox_p <- suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  out
}

# DT column renderer: fixed decimals or 3 significant digits; missing values shown as "-"
dt_num <- function(digits = 3, signif = FALSE) {
  fmt <- if (signif)
    "(Math.abs(x) >= 1e-3 && Math.abs(x) < 1e4) ? String(Number(x.toPrecision(3))) : x.toExponential(2)"
  else sprintf("x.toFixed(%d)", digits)
  JS(sprintf(paste0(
    "function(d, type) { if (type !== 'display') return d;",
    " if (d === null || d === '' || isNaN(Number(d))) return '-';",
    " var x = Number(d); return %s; }"), fmt))
}

# ---- Correlation helpers ----------------------------------------------------
# Pearson and Spearman correlation of each feature with `target`, using the bins
# where both values are finite. Features with < 3 such bins or no variance get NA.
cor_with_target <- function(d, features, target) {
  y <- as.numeric(d[[target]])
  res <- vapply(features, function(f) {
    x  <- as.numeric(d[[f]])
    ok <- is.finite(x) & is.finite(y)
    n  <- sum(ok)
    if (n < 3 || !isTRUE(var(x[ok]) > 0) || !isTRUE(var(y[ok]) > 0)) return(c(NA, NA, n))
    c(cor(x[ok], y[ok]), cor(x[ok], y[ok], method = "spearman"), n)
  }, numeric(3))
  out <- data.frame(Feature = features, Pearson = res[1, ], Spearman = res[2, ],
                    n = as.integer(res[3, ]), stringsAsFactors = FALSE, row.names = NULL)
  out[order(-abs(out$Spearman), -abs(out$Pearson), na.last = TRUE), , drop = FALSE]
}

# Numeric matrix of `vars`; drops variables with < 3 finite values or no variance
numeric_matrix <- function(d, vars) {
  M <- vapply(vars, function(v) as.numeric(d[[v]]), numeric(nrow(d)))
  if (!is.matrix(M)) M <- matrix(M, nrow = 1, dimnames = list(NULL, vars))
  M[!is.finite(M)] <- NA
  keep <- apply(M, 2, function(x) sum(!is.na(x)) >= 3 && isTRUE(var(x, na.rm = TRUE) > 0))
  list(M = M[, keep, drop = FALSE], dropped = vars[!keep])
}

# Correlation matrix on pairwise-complete bins; same result as
# cor(M, method = method, use = "pairwise.complete.obs"). For Spearman, base R
# re-ranks both columns for every pair, which is slow with many features over
# all chromosomes; here a column's ranks are reused whenever the pair loses
# none of that column's bins.
cor_matrix <- function(M, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (method == "pearson")
    return(suppressWarnings(cor(M, use = "pairwise.complete.obs")))
  p  <- ncol(M)
  ok <- !is.na(M)
  rk <- apply(M, 2, rank, na.last = "keep")
  if (all(ok)) return(suppressWarnings(cor(rk)))
  R  <- diag(1, p); dimnames(R) <- list(colnames(M), colnames(M))
  ranks_on <- function(k, both)    # reuse the column's ranks when no bins are lost
    if (identical(both, ok[, k])) rk[both, k] else rank(M[both, k])
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    both <- ok[, i] & ok[, j]
    R[i, j] <- R[j, i] <- if (sum(both) < 3) NA_real_
    else suppressWarnings(cor(ranks_on(i, both), ranks_on(j, both)))
  }
  R
}

# Variable order for the matrix: hierarchical clustering on 1 - r
# (missing correlations treated as 0), or the given order otherwise.
cluster_order <- function(R) {
  if (ncol(R) < 3) return(seq_len(ncol(R)))
  R0 <- R; R0[is.na(R0)] <- 0
  hclust(as.dist(1 - R0), method = "average")$order
}


# ---- Scatterplot matrix ------------------------------------------------------
# Lower panels: scatter + lowess (red); diagonal: density + variable name;
# upper panels: correlation coefficient with significance stars (red).
# Coefficients and p-values use every bin with both values (pairwise); only the
# scatter plots are drawn from a random subset of bins, to keep it fast.
PAIRS_MAX_POINTS <- 5000

cor_stars <- function(p) {
  if (!is.finite(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**"
  else if (p < 0.05) "*" else ""
}

# Five evenly spaced tick labels across the variable's range, rotated so they
# never overlap; trailing zeros are dropped (e.g. 0.25 rather than 0.2500).
pairs_axis <- function(side, r) {
  at <- seq(r[1], r[2], length.out = 5)
  axis(side, at = at, labels = formatC(signif(at, 3), format = "fg", digits = 3),
       las = 2, cex.axis = 0.85, col = "grey40", col.ticks = "grey40")
}

# Colours of the two groups when bins are split by a feature threshold
SPLIT_COLORS <- c(present = "#009E73", absent = "#E69F00")

# Correlation p-values from r and the pairwise n (t-approximation, as
# cor.test() uses for Pearson and non-exact Spearman)
cor_pvals <- function(R, M) {
  N <- crossprod(!is.na(M))
  tt <- R * sqrt((N - 2) / pmax(1 - R^2, 1e-12))
  2 * pt(-abs(tt), df = pmax(N - 2, 1))
}

# M: bins x variables; R: correlation matrix of M. Optional g: logical vector
# (TRUE = feature present, FALSE = absent, NA = unknown) splitting the bins
# into two groups that get their own colour, lowess curve, density and
# coefficient; labs = the two group labels (present, absent).
draw_pairs <- function(M, R, method = "spearman", g = NULL, labs = NULL) {
  p <- ncol(M); v <- colnames(M)
  split <- !is.null(g)
  if (split) {
    grp <- list(present = which(g %in% TRUE), absent = which(g %in% FALSE))
    Rg <- lapply(grp, function(ix) if (length(ix) >= 3)
      suppressWarnings(cor_matrix(M[ix, , drop = FALSE], method)) else R * NA)
    Pg <- Map(function(Rk, ix) cor_pvals(Rk, M[ix, , drop = FALSE]), Rg, grp)
  } else P <- cor_pvals(R, M)
  set.seed(1)
  rows <- if (nrow(M) > PAIRS_MAX_POINTS) sort(sample(nrow(M), PAIRS_MAX_POINTS))
          else seq_len(nrow(M))
  S <- M[rows, , drop = FALSE]
  gs <- if (split) g[rows] else NULL
  lab <- vapply(gsub("_", " ", v),
                function(z) paste(strwrap(z, width = 12), collapse = "\n"), character(1))
  rng <- lapply(seq_len(p), function(j) range(M[, j], na.rm = TRUE))
  cex_num <- if (p > 15) 1.1 else if (p > 10) 1.5 else 1.9
  if (split) cex_num <- cex_num * 0.75
  cex_lab <- if (p > 15) 0.75 else 0.95
  pt_col  <- if (split) ifelse(gs %in% TRUE, adjustcolor(SPLIT_COLORS[["present"]], 0.5),
                         ifelse(gs %in% FALSE, adjustcolor(SPLIT_COLORS[["absent"]], 0.5), NA))
             else adjustcolor("black", alpha.f = 0.45)
  smooth <- function(x, y, col) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) > 10 && length(unique(x[ok])) > 2) lines(lowess(x[ok], y[ok]), col = col, lwd = 2)
  }
  dens <- function(x) { x <- x[is.finite(x)]; if (length(unique(x)) > 1) density(x) else NULL }

  op <- par(mfrow = c(p, p), mar = rep(0.25, 4),
            oma = c(4.5, 4.5, if (split) 2.2 else 0.5, 0.5),
            mgp = c(2, 0.45, 0), tcl = -0.25, col.axis = "grey15")
  on.exit(par(op))
  for (i in seq_len(p)) for (j in seq_len(p)) {
    if (i == j) {                                     # diagonal: density + name
      dl <- if (split) lapply(grp, function(ix) dens(M[ix, j])) else list(all = dens(M[, j]))
      ymax <- max(c(1e-12, unlist(lapply(dl, function(dd) if (!is.null(dd)) max(dd$y)))))
      plot(NA, xlim = rng[[j]], ylim = c(0, ymax * 1.35),
           axes = FALSE, xlab = "", ylab = "", xaxs = "i")
      for (k in names(dl)) {
        dd <- dl[[k]]; if (is.null(dd)) next
        col <- if (split) SPLIT_COLORS[[k]] else "black"
        polygon(c(dd$x, rev(dd$x)), c(dd$y, rep(0, length(dd$y))),
                col = if (split) adjustcolor(col, 0.2) else "grey90", border = NA)
        lines(dd, col = col)
      }
      text(mean(rng[[j]]), par("usr")[4] * 0.72, lab[j], col = "darkblue",
           font = 2, cex = cex_lab)
    } else if (i > j) {                               # lower: scatter + lowess
      x <- S[, j]; y <- S[, i]
      plot(x, y, xlim = rng[[j]], ylim = rng[[i]], pch = 16, cex = 0.55, col = pt_col,
           axes = FALSE, xlab = "", ylab = "")
      if (split) {
        smooth(x[gs %in% TRUE],  y[gs %in% TRUE],  SPLIT_COLORS[["present"]])
        smooth(x[gs %in% FALSE], y[gs %in% FALSE], SPLIT_COLORS[["absent"]])
      } else smooth(x, y, "red")
    } else {                                          # upper: coefficient + stars
      plot.new(); plot.window(c(0, 1), c(0, 1))
      fmt <- function(r) if (is.finite(r)) sprintf("%.2f", r) else "NA"
      if (split) {
        for (k in c("present", "absent")) {
          yy <- if (k == "present") 0.68 else 0.32
          text(0.5, yy, paste0(fmt(Rg[[k]][i, j]), cor_stars(Pg[[k]][i, j])),
               col = SPLIT_COLORS[[k]], font = 2, cex = cex_num)
        }
      } else {
        text(0.5, 0.5, fmt(R[i, j]), font = 2, cex = cex_num)
        text(0.97, 0.9, cor_stars(P[i, j]), col = "red", font = 2, adj = c(1, 0.5),
             cex = cex_num * 0.75)
      }
    }
    box(col = "grey20")
    if (i == p)          pairs_axis(1, rng[[j]])
    if (j == 1 && i > 1) pairs_axis(2, rng[[i]])
  }
  if (split) {                                        # legend across the top
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("top", legend = sprintf("%s (n = %s)", labs,
                                   format(lengths(grp), big.mark = ",")),
           col = SPLIT_COLORS, pch = 16, lwd = 2, horiz = TRUE, bty = "n",
           cex = 0.95, text.col = SPLIT_COLORS)
  }
}

# ---- Slide-out documentation panel ------------------------------------------
doc_css <- "
.doc-drawer{position:fixed;top:0;right:0;height:100%;width:520px;max-width:92vw;
  background:#fff;box-shadow:-4px 0 18px rgba(0,0,0,.25);z-index:1040;
  transform:translateX(100%);transition:transform .35s ease;}
.doc-drawer.open{transform:translateX(0);}
.doc-tab{position:absolute;left:-38px;top:35%;width:38px;padding:14px 0;
  writing-mode:vertical-rl;text-align:center;cursor:pointer;user-select:none;
  background:#2c3e50;color:#fff;border-radius:8px 0 0 8px;font-weight:600;
  letter-spacing:.05em;box-shadow:-3px 2px 8px rgba(0,0,0,.2);}
.doc-tab:hover{background:#1a252f;}
.doc-body{height:100%;overflow-y:auto;padding:20px 24px 60px 24px;font-size:13px;}
.doc-body h3{margin-top:0;}
.doc-body h4{margin-top:22px;border-bottom:1px solid #ddd;padding-bottom:4px;}
.doc-close{float:right;cursor:pointer;font-size:22px;line-height:1;color:#666;}
"

doc_drawer <- function() {
  tagList(
    tags$style(HTML(doc_css)),
    div(id = "doc-drawer", class = "doc-drawer open",
        div(class = "doc-tab",
            onclick = "document.getElementById('doc-drawer').classList.toggle('open')",
            "Documentation"),
        div(class = "doc-body",
            span(class = "doc-close", HTML("&times;"),
                 onclick = "document.getElementById('doc-drawer').classList.remove('open')"),
            h3("Feature documentation"),
            
            h4("Genome-specific feature pairs"),
            p("For each donor parent, features are provided for each parental genome involved in",
              "the cross (P1: MbTMV; P2: Pen, Che or MT). Coordinates were converted between the",
              "two genomes as described in \"Coordinate conversion\" below. For example, \"GA_P1\"",
              "holds the genome annotations of the recurrent parent (RP) genome, MbTMV, whereas",
              "\"GA_P2\" holds the annotations of the donor parent (DP) genome projected onto the RP",
              "coordinates. This applies to every feature except SVs, which depend implicitly on",
              "the two crossing parents."),
            
            h4("Coordinate conversion"),
            p("Coordinates were converted with a custom script (syri_liftover.py) that combines the",
              "minimap2 BAM alignment between the two parental genomes with the SyRI annotation of",
              "that same alignment. The start and end of every interval are converted separately,",
              "and only syntenic (SYN/SYNAL) and inverted (INV/INVAL) conversions are kept:"),
            tags$ul(
              tags$li(strong("Exact"), " – the position lies in an aligned block of a syntenic",
                      " or inverted alignment and is converted base-to-base using the alignment",
                      " (CIGAR). Alignments that SyRI did not report but that lie completely inside",
                      " a SYN or INV region on both genomes are also used."),
              tags$li(strong("Interpolated"), " – the position lies inside a SYN or INV region",
                      " but in a stretch without alignment (e.g. a highly divergent region or a",
                      " deletion). Its coordinate is interpolated linearly between the closest",
                      " aligned bases on either side within the region."),
              tags$li(strong("Not converted"), " – positions outside every SYN and INV region",
                      " (translocations, duplications, unaligned regions) get no coordinate, and",
                      " the interval is dropped. The exception is centromeric regions (see",
                      " \"Syntenic fallback\").")),
            p("Intervals in inversions are reported with start < end in the target genome.",
              "When several syntenic or inverted conversions exist for one position, the one on",
              "the same chromosome closest to the original position is kept. Positions within",
              "50 bp of a chromosome end are mapped to the corresponding end of the homologous",
              "chromosome. Every candidate conversion (including translocated and duplicated",
              "alignments that were not used) is recorded in a log file next to each feature",
              "file (*_conversion_hits.tsv.gz), together with its SyRI annotation and the method",
              "used."),
            
            h4("Syntenic fallback"),
            p("For centromeric regions only, a boundary that falls outside every SYN and INV",
              "region is moved to the closest syntenic or inverted aligned base: upstream for the",
              "interval start and downstream for the interval end. That base is then converted",
              "exactly. The rationale is that every genome must contain at least one centromeric",
              "region, and expanding the interval because of the missing direct mapping is not",
              "expected to have a large impact, since tomato centromeres are surrounded by long",
              "stretches of heterochromatin, which also inhibits CO formation."),
            
            h4("Structural variant features"),
            p("Structural variant features were extracted from SyRI (v1.8.2) annotations of",
              "genome-wide minimap2 (v2.31-r1302) alignments between the reference genomes of the",
              "parents in each cross. Sequence divergence was accommodated per species (-x asm5 for",
              "MicroTom, -x asm20 for Cheesmaniae and Pennellii). SV features comprise three",
              "high-order features (SYN, INV and NLO) and five low-order features (SYNAL, HDR, TDM,",
              "INV and NLO). High-order features are categorized by the expected outcome of",
              "recombination:"),
            tags$ul(
              tags$li(strong("Syntenic regions (SYN)"), " align properly during synapsis because",
                      " they are locally orthologous between the homologous chromosomes, which",
                      " enables recombination. They come from SyRI annotation SYN (column 11 of the",
                      " SyRI .out file) and may contain syntenic alignments (SYNAL), highly",
                      " divergent regions (HDR) or tandem repeats (TDM)."),
              tags$li(strong("Inversions (INV)"), " align by forming an inverted loop, so a gamete",
                      " resulting from any crossover within the region is non-viable. Inversions",
                      " correspond to the value INV in column 11 of the SyRI .out file."),
              tags$li(strong("Non-locally orthologous regions (NLO)"), " lack a 1-to-1 orthologous",
                      " sequence at the same location on the homologous pair. For example,",
                      " translocated copies have an ortholog elsewhere in the homologous genome, but",
                      " not at the same location, so they cannot align properly during synapsis and",
                      " are therefore non-recombinogenic. NLO regions contain the annotations",
                      " translocations (TRANS), inverted translocations (INVTR), copy loss (CPL),",
                      " duplications (DUP) with copy loss in the reference, and inverted duplications",
                      " (INVDP) with copy loss in the reference (MbTMV is the reference). SVs with a",
                      " copy gain in the query (Che, Pen, MT) do not interrupt synteny in the MbTMV",
                      " coordinate system and were discarded. In these cases, the MbTMV sequences",
                      " flanking the insertion point are neighbors and can align with their homologs",
                      " despite the intervening SV, whereas only the insertion in the homolog cannot",
                      " recombine, because no homologous sequence exists in MbTMV.")),
            
            h4("Feature categories"),
            tags$ul(
              tags$li(strong("Genome (GA)"), " – genome annotations, as the fraction of the",
                      " window covered."),
              tags$li(strong("Transposable elements (TE)"), " – transposable element",
                      " annotations, as the fraction of the window covered."),
              tags$li(strong("Chromatin"), " – chromatin compartments, given as Hi-C PCA",
                      " eigenvectors from -1 to 1 (positive = active chromatin)."),
              tags$li(strong("Repeats"), " – repeat annotation score; high values indicate",
                      " repetitive regions."),
              tags$li(strong("Centromeres (CEN)"), " – centromeric regions, as the fraction of",
                      " the window covered."),
              tags$li(strong("Relative chromosome position (POS)"), " – position of the window",
                      " along its chromosome arm, from 0 (telomere) to 1 (centromere).")),
            p("Structural variant (SV) features are described above."),
            
            h4("Interpretation of the feature values"),
            p("The genome annotations (GA) and transposable elements (TE) were taken from the",
              "annotations presented by Rengs et al., 2026. The centromeric regions were provided",
              "by Elias. All these features (GA*, TE*, CEN*, CHR*) and their low-order features are",
              "quantified as the fraction of the genomic window covered by the corresponding",
              "landmark. For example, CEN_P1 = 0.80 means that 80% of the window in the recurrent",
              "parent genome is covered by centromeric regions."),
            p("The chromatin compartments were also provided by Elias. They are quantified as the",
              "eigenvectors of a PCA of the Hi-C distance matrix. Signs were assigned based on the",
              "correlation with gene content, so values span -1 to 1, with positive values",
              "reflecting active chromatin. Repeat annotations (Rengs et al., 2026) are quantified",
              "as a score for the likelihood of a repetitive region, with high values meaning",
              "repeats."),
            p("The relative chromosome position (POS_P1, POS_P2) is computed from the window",
              "midpoint and the centromere defined by CEN_P1 or CEN_P2, separately for each",
              "chromosome and cross. On the left arm, it is the distance from the start of the",
              "chromosome divided by the distance from there to the midpoint of the first window",
              "with CEN = 1; on the right arm, the same is measured from the end of the chromosome",
              "back to the midpoint of the last window with CEN = 1. All windows between the first",
              "and the last CEN = 1 window are set to 1, so a patchy centromere counts as one block.",
              "Because the telomere is taken as the start of the first window, the first window",
              "never reaches exactly 0. Chromosomes without any window with CEN = 1 get NA."),
            
            h4("Notes"),
            tags$ul(
              tags$li("When collecting data, I assumed that intervals were 0-based half-open for BED",
                      " files and 1-based closed for GFF/GFF3 files."),
              tags$li("Features are only converted through syntenic and inverted regions, so",
                      " features of the donor genome that sit in translocations or duplications are",
                      " absent from the P2 features. The SV features (INV, NLO) still describe these",
                      " regions, so the ML models can learn that they are strong predictors of no",
                      " recombination. In fact, the INV/NLO signal appears very clear in coldspots."),
              tags$li("This app and its content are at an early stage. Any comment or idea for",
                      " improvement is very welcome."))
        )
    )
  )
}
# ---- UI ---------------------------------------------------------------------
corr_css <- "
.stat-panel table.dataTable{font-size:12px;}
.stat-panel table.dataTable td, .stat-panel table.dataTable th{padding:3px 6px;white-space:nowrap;}
.stat-panel .dataTables_filter{float:none;text-align:left;}
.stat-panel .dataTables_filter input{width:100%;margin-left:0;}
"

corr_panel <- wellPanel(
  id = "corr_panel", class = "stat-panel",
  tags$style(HTML(corr_css)),
  strong("Correlation with recombination rate"),
  selectInput("corr_target", "Recombination rate", choices = NULL),
  radioButtons("corr_scope", "Bins used",
               choices = c("Current chromosome" = "chr", "All chromosomes" = "all"),
               inline = TRUE),
  DTOutput("corr_table"),
  uiOutput("corr_note"),
  actionButton("corr_open", "Open correlation matrix", icon = icon("table"),
               class = "btn-primary btn-sm", style = "margin-top:6px;")
)

poi_panel <- wellPanel(
  id = "poi_panel", class = "stat-panel",
  strong("Feature vs. CO state (Poisson test)"),
  selectizeInput("poi_feature", "Feature", choices = NULL),
  sliderInput("poi_thr", "Feature present if value >", min = 0, max = 1, value = 0.5,
              step = 0.01),
  helpText("Default: 0.50 for inversions and chromatin state, median for other features."),
  selectInput("poi_state", "CO state",
              choices = c("No CO in bin" = "noco", "\u2265 1 CO in bin" = "co",
                          "Coldspot" = "cold", "Hotspot" = "hot")),
  radioButtons("poi_scope", "Chromosomes",
               choices = c("Chromosome-wise" = "chr", "All chromosomes" = "all"),
               inline = TRUE),
  numericInput("poi_ntests", "Bonferroni: number of features screened", value = 1,
               min = 1, step = 1),
  DTOutput("poi_table"),
  div(id = "poi_warning",
      style = paste("margin:8px 0;padding:8px 10px;background:#fff4e5;",
                    "border-left:4px solid #e69f00;border-radius:3px;color:#5a3d00;"),
      strong("Keep the Poisson test for its effect size (O/E), but don't rely on its",
             "p-values: the many zero bins and the hotspots make the variance much larger",
             "than the Poisson model assumes. Use the quasi-Poisson p-value instead.")),
  uiOutput("poi_note")
)

main_ui <- fluidRow(
  # Left column: controls, with the correlation panel below them
  column(
    width = 3,
    wellPanel(
      fileInput("files", "Load BED files (one or more resolutions)", multiple = TRUE,
                accept = c(".bed", ".txt", ".tsv")),
      helpText("Files are read in your browser and never leave your computer."),
      selectInput("chr", "Chromosome (Chr_ID)", choices = NULL),
      selectInput("donor", "Donor parent", choices = NULL),
      selectInput("male", "Male (meiosis)", choices = NULL),
      helpText("Male values: 0 = female meiosis, 1 = male meiosis, NA = both sexes combined."),
      selectInput("dataset", "Dataset", choices = NULL),
      selectizeInput("ys", "Variables to plot (y)", choices = NULL, multiple = TRUE),
      div(style = "margin-bottom:10px;",
          actionButton("ys_all", "Select all", class = "btn-sm"),
          actionButton("ys_none", "Clear variables", class = "btn-sm")),
      actionButton("dist_open", "Show distributions & correlations",
                   icon = icon("chart-column"), class = "btn-primary btn-sm",
                   style = "margin-bottom:10px;"),
      radioButtons("layout", "Layout",
                   choices = c("Separate panels (free y)" = "panels",
                               "Overlay (each variable scaled 0-1)" = "overlay")),
      hr(),
      strong("Hotspots / coldspots"),
      checkboxGroupInput("show_types", NULL, choices = c("Hotspot", "Coldspot"),
                         selected = c("Hotspot", "Coldspot"), inline = TRUE),
      numericInput("min_len", "Minimum region length (bp)", 500000, min = 1, step = 50000),
      numericInput("hot_p", "Hotspot: mean Poisson_p <", 0.1, min = 0, max = 1, step = 0.01),
      numericInput("cold_p", "Coldspot: mean Poisson_p >=", 1, min = 0, max = 1, step = 0.01),
      helpText("Regions are runs of consecutive bins spanning >= the minimum length whose",
               "mean Poisson_p meets the threshold; overlapping runs are merged.",
               "Hotspots are shaded red, coldspots blue."),
      hr(),
      helpText("Click points to select them (click again to deselect). Use the box/lasso",
               "tool in the plot toolbar to select many."),
      actionButton("clear", "Clear selection")
    ),
    poi_panel,
    corr_panel
  ),
  column(
    width = 9,
    tabsetPanel(
      id = "main_tabs",
      tabPanel(
        "Line plots", value = "lines",
        h4("Crossover statistics for the current selection"),
        tableOutput("stats"),
        uiOutput("stats_note"),
        hr(),
        plotlyOutput("plot", height = "auto"),
        hr(),
        textOutput("sel_info"),
        DTOutput("table")
      ),
      tabPanel(
        "Distributions & correlations", value = "dist",
        h4("Distributions and pairwise correlations of the plotted variables"),
        fluidRow(
          column(4, radioButtons("dv_scope", "Bins used",
                                 c("Current chromosome" = "chr", "All chromosomes" = "all"),
                                 inline = TRUE)),
          column(4, radioButtons("dv_method", "Correlation",
                                 c("Spearman" = "spearman", "Pearson" = "pearson"),
                                 inline = TRUE)),
          column(4, radioButtons("dv_order", "Correlation order",
                                 c("As selected" = "original",
                                   "Hierarchical clustering" = "cluster"),
                                 inline = TRUE))
        ),
        fluidRow(
          column(4, selectizeInput("dv_split", "Split bins by feature", choices = NULL)),
          column(4, conditionalPanel(
            "input.dv_split && input.dv_split !== 'none'",
            sliderInput("dv_split_thr", "Feature present if value >", min = 0, max = 1,
                        value = 0.5, step = 0.01, width = "100%")))
        ),
        uiOutput("dv_note"),
        hr(),
        h4("Distributions"),
        plotOutput("dv_box", height = "auto"),
        hr(),
        h4("Pairwise correlations \u2013 all bins"),
        uiOutput("dv_pairs_note"),
        plotOutput("dv_pairs", height = "auto"),
        hr(),
        h4("Pairwise correlations \u2013 recombinogenic bins only (\u2265 1 CO)"),
        uiOutput("dv_pairs_rec_note"),
        plotOutput("dv_pairs_rec", height = "auto")
      )
    )
  )
)

ui <- fluidPage(
  titlePanel("The recombination landscape of tomato crosses – An exploratory analysis"),
  tags$button(
    type = "button", class = "btn btn-default btn-sm", style = "margin-bottom:10px;",
    onclick = "document.getElementById('doc-drawer').classList.toggle('open')",
    "\U0001F4D6 Show / hide documentation"),
  main_ui,
  doc_drawer()   # slide-out documentation
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  
  raw <- reactive({
    validate(need(!is.null(input$files),
                  "Load one or more BED files (left panel) to start."))
    paths <- input$files$datapath
    d <- withProgress(message = "Reading files...", value = 0.3,
                      combine_tables(lapply(paths, read_bed)))
    d[, Interval_midpoint := (Interval_start + Interval_end) / 2]
    d
  })
  
  feature_cols <- reactive({
    d <- raw()
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    list(targets  = intersect(TARGET_COLS, names(d)),
         features = setdiff(num, c(META_COLS, TARGET_COLS)))
  })
  
  # Populate controls when data are loaded
  observeEvent(raw(), {
    d <- raw()
    updateSelectInput(session, "chr", choices = sort(unique(d$Chr_ID)))
    updateSelectInput(session, "donor", choices = sort(unique(d$Donor_parent)))
    updateSelectInput(session, "dataset", choices = sort(unique(d$Dataset)))
    lab <- c("NA" = "Both sexes (NA)", "0" = "Female (0)", "1" = "Male (1)")
    present <- unique(ifelse(is.na(d$Male), "NA", as.character(d$Male)))
    present <- intersect(names(lab), present)
    updateSelectInput(session, "male", choices = setNames(present, lab[present]))
    fc <- feature_cols()
    cats <- factor(vapply(fc$features, feature_category, character(1), USE.NAMES = FALSE),
                   levels = CATEGORY_LEVELS)
    ch <- c(list(Targets = fc$targets), split(fc$features, cats))
    ch <- lapply(ch[lengths(ch) > 0], as.list)   # as.list: keeps one-item groups as groups
    sel <- intersect(DEFAULT_PLOT_VARS, c(fc$targets, fc$features))
    if (!length(sel)) sel <- c(head(fc$targets, 1), head(fc$features, 2))
    updateSelectizeInput(session, "ys", choices = ch, selected = sel)
    fch <- split(fc$features, cats)
    fch <- lapply(fch[lengths(fch) > 0], as.list)
    curs <- isolate(input$dv_split)
    updateSelectizeInput(session, "dv_split", choices = c(list("No split" = "none"), fch),
                         selected = if (isTRUE(curs %in% fc$features)) curs else "none")
    inv <- Filter(function(f) has_token(f, "INV"), fc$features)
    curf <- isolate(input$poi_feature)
    updateSelectizeInput(session, "poi_feature", choices = fch,
                         selected = if (isTRUE(curf %in% fc$features)) curf
                         else c(inv, fc$features)[1])
    cur <- isolate(input$corr_target)
    updateSelectInput(session, "corr_target", choices = fc$targets,
                      selected = if (isTRUE(cur %in% fc$targets)) cur else head(fc$targets, 1))
  })
  
  observeEvent(input$ys_all, {
    fc <- feature_cols()
    updateSelectizeInput(session, "ys", selected = c(fc$targets, fc$features))
  })
  observeEvent(input$ys_none, {
    updateSelectizeInput(session, "ys", selected = character(0))
  })
  
  selected <- reactiveVal(integer(0))
  observeEvent(list(input$chr, input$donor, input$male, input$dataset),
               selected(integer(0)), ignoreInit = TRUE)
  observeEvent(input$clear, selected(integer(0)))
  
  # Filtered data for the current Chr / donor / sex / dataset + spot calling
  dat <- reactive({
    req(input$chr, input$donor, input$male, input$dataset)
    d <- subset_rows(raw(), input$donor, input$male, input$dataset)
    d <- d[Chr_ID == input$chr]
    validate(need(nrow(d) > 0, "No rows for this combination of selections."))
    setorder(d, Interval_start)
    d[, .id := .I]
    sp <- call_spots(d, input$min_len, input$hot_p, input$cold_p)
    d[, Spot_class := sp$class]
    list(d = d, regions = sp$regions)
  })
  
  # ---- CO statistics (chromosome x donor/species x sex) ----
  # Bins overlap (sliding window), so summing Number_of_COs over all bins would
  # count most COs twice. We therefore sum over a non-overlapping tiling:
  # greedily keep each bin that starts at or after the end of the last kept bin.
  co_stats <- reactive({
    d <- dat()$d
    validate(need("Number_of_COs" %in% names(d) && "Number_of_hybrids" %in% names(d),
                  "Number_of_COs / Number_of_hybrids columns not found."))
    tiled <- d[tile_keep(Interval_start, Interval_end)]
    gam <- unique(na.omit(d$Number_of_hybrids))
    n_gam <- max(gam)
    total <- sum(tiled$Number_of_COs, na.rm = TRUE)
    list(chr = input$chr, donor = input$donor, male = input$male,
         gametes = n_gam, gametes_vary = length(gam) > 1,
         total = total, per_gamete = total / n_gam,
         n_bins = nrow(tiled), n_all = nrow(d),
         chr_len_mb = max(d$Interval_end) / 1e6)
  })
  
  output$stats <- renderTable({
    s <- co_stats()
    data.frame(
      Chromosome = s$chr,
      `Donor parent` = s$donor,
      Meiosis = c("0" = "Female", "1" = "Male", "NA" = "Both sexes")[[s$male]],
      `Gametes (hybrids)` = format(s$gametes, big.mark = ","),
      `Total COs` = format(s$total, big.mark = ","),
      `COs per gamete` = sprintf("%.3f", s$per_gamete),
      `COs per gamete per Mb` = sprintf("%.4f", s$per_gamete / s$chr_len_mb),
      `cM/Mbp` = sprintf("%.3f", 100 * s$per_gamete / s$chr_len_mb),
      check.names = FALSE
    )
  }, align = "c")
  
  output$stats_note <- renderUI({
    s <- co_stats()
    helpText(sprintf(
      paste0("Each hybrid = one gamete. Total COs are summed over %d non-overlapping bins ",
             "(of %d sliding-window bins) to avoid double counting; a partial bin at the ",
             "chromosome end may be omitted, so the total is approximate. cM/Mbp = COs per gamete per Mb x 100 ",
             "(1 CO per gamete = 1 Morgan = 100 cM).%s"),
      s$n_bins, s$n_all,
      if (s$gametes_vary) " Number_of_hybrids varies between bins; the maximum is used." else ""))
  })
  
  # ---- Selection events ----
  # Before the plot is first drawn plotly warns that the events are not registered
  # yet; harmless, so silenced (the plot registers them with event_register())
  click <- reactive(suppressWarnings(
    event_data("plotly_click", source = "main", priority = "event")))
  observeEvent(click(), {
    ids <- suppressWarnings(as.integer(unlist(click()$key)))
    ids <- ids[!is.na(ids)]
    if (!length(ids)) return()
    cur <- selected()
    selected(if (all(ids %in% cur)) setdiff(cur, ids) else union(cur, ids))
  })
  brush <- reactive(suppressWarnings(event_data("plotly_selected", source = "main")))
  observeEvent(brush(), {
    ids <- suppressWarnings(as.integer(unlist(brush()$key)))
    ids <- ids[!is.na(ids)]
    if (length(ids)) selected(union(selected(), ids))
  })
  
  # ---- Plot ----
  output$plot <- renderPlotly({
    x <- dat(); d <- x$d
    ys <- input$ys
    validate(need(length(ys) > 0, "Select at least one variable to plot."))
    overlay <- input$layout == "overlay"
    
    # Cross label: Background_parent x Donor_parent
    bg <- if ("Background_parent" %in% names(d))
      paste(sort(unique(na.omit(d$Background_parent))), collapse = "/") else ""
    cross <- if (nzchar(bg)) paste(bg, "x", input$donor) else input$donor
    
    long <- melt(d, id.vars = c(".id", "Interval_midpoint", intersect("Number_of_COs", names(d))),
                 measure.vars = ys,
                 variable.name = "variable", value.name = "value",
                 variable.factor = FALSE)
    long[, variable := factor(variable, levels = ys)]
    long[, mid := Interval_midpoint / 1e6]
    if (overlay) long[, plot_value := rescale01(value), by = variable] else long[, plot_value := value]
    long[, hover := paste0(variable, "<br>", round(mid, 3), " Mb<br>value: ", signif(value, 4))]
    if ("Number_of_COs" %in% names(long))   # CO count in the tooltip of recombination-rate traces
      long[variable %in% TARGET_COLS, hover := paste0(hover, "<br>COs: ", Number_of_COs)]
    
    p <- ggplot()
    
    # Hot/cold spot shading
    reg <- copy(x$regions)[type %in% input$show_types]
    if (nrow(reg) > 0) {
      rng <- if (overlay) data.table(ymin = 0, ymax = 1)
      else long[, .(ymin = min(plot_value, na.rm = TRUE),
                    ymax = max(plot_value, na.rm = TRUE)), by = variable]
      rng[ymin == ymax, `:=`(ymin = ymin - 0.5, ymax = ymax + 0.5)]
      rng[, k := 1]; reg[, k := 1]
      rects <- merge(reg, rng, by = "k", allow.cartesian = TRUE)
      rects[, `:=`(xmin = start / 1e6, xmax = end / 1e6)]
      p <- p + geom_rect(data = rects,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
                         alpha = 0.2, inherit.aes = FALSE)
    }
    
    # key / text are plotly-only aesthetics: ggplot2 warns about them, ggplotly uses them
    p <- p +
      geom_line(data = long, aes(mid, plot_value, colour = variable, group = variable)) +
      suppressWarnings(
        geom_point(data = long, aes(mid, plot_value, colour = variable, key = .id, text = hover),
                   size = 0.5))
    
    sel <- selected()
    if (length(sel)) {
      p <- p + suppressWarnings(
        geom_point(data = long[.id %in% sel], aes(mid, plot_value, key = .id),
                   shape = 21, size = 4, stroke = 1.2, colour = "black", fill = NA,
                   inherit.aes = FALSE))
    }
    
    p <- p +
      scale_fill_manual(values = SPOT_COLORS, name = "Region") +
      labs(x = "Interval_midpoint (Mb)",
           y = if (overlay) "Value (min-max scaled per variable)" else NULL,
           colour = NULL,
           title = sprintf("%s | %s | %s", input$chr, cross,
                           c("0" = "female meiosis", "1" = "male meiosis",
                             "NA" = "both sexes")[[input$male]])) +
      theme_minimal(base_size = 12)
    if (!overlay) p <- p + facet_wrap(~variable, ncol = 1, scales = "free_y")
    
    h <- if (overlay) 450 else max(350, 210 * length(ys))
    pl <- suppressWarnings(plotly_build(
      ggplotly(p, tooltip = "text", source = "main", height = h)))
    pl <- merge_spot_legend(pl)
    pl %>% layout(uirevision = "keep", dragmode = "zoom") %>%
      event_register("plotly_click") %>%
      event_register("plotly_selected")
  })
  
  # ---- Correlation analysis ----
  # Bins: the current chromosome or all chromosomes, for the selected donor /
  # sex / dataset. Every sliding-window bin is used.
  corr_data <- reactive({
    req(input$donor, input$male, input$dataset, input$corr_scope)
    d <- subset_rows(raw(), input$donor, input$male, input$dataset)
    if (input$corr_scope == "chr") {
      req(input$chr)
      d <- d[Chr_ID == input$chr]
    }
    validate(need(nrow(d) >= 3, "Too few bins for correlation."))
    d
  })
  
  corr_scope_label <- reactive({
    sex <- c("0" = "female meiosis", "1" = "male meiosis", "NA" = "both sexes")[[input$male]]
    where <- if (input$corr_scope == "chr") input$chr else "all chromosomes"
    paste(where, input$donor, sex, input$dataset, sep = " | ")
  })
  
  corr_table <- reactive({
    req(input$corr_target)
    d <- corr_data(); fc <- feature_cols()
    validate(need(input$corr_target %in% names(d), "Recombination-rate column not found."),
             need(length(fc$features) > 0, "No feature columns found."))
    cor_with_target(d, fc$features, input$corr_target)
  })
  
  output$corr_table <- renderDT({
    ct <- corr_table()
    cuts <- round(seq(-0.95, 0.95, by = 0.1), 2)            # middle bin (-0.05, 0.05] is neutral
    cols <- colorRampPalette(CORR_TABLE_POLES)(length(cuts) + 1)
    datatable(ct, rownames = FALSE, selection = "none", class = "compact hover",
              colnames = c("Feature", "Pearson r", "Spearman \u03c1", "n"),
              options = list(dom = "ft", paging = FALSE, scrollY = "320px",
                             scrollX = TRUE, scrollCollapse = TRUE, order = list(),
                             columnDefs = list(list(targets = 1:2, render = dt_num(3)),
                                               list(targets = 3, render = dt_num(0))),
                             language = list(search = "",
                                             searchPlaceholder = "Filter features"))) %>%
      formatStyle(c("Pearson", "Spearman"), backgroundColor = styleInterval(cuts, cols))
  })
  
  output$corr_note <- renderUI({
    n <- nrow(corr_data()); req(corr_table())
    helpText(sprintf(paste0(
      "%s bins (%s). Sorted by |Spearman \u03c1|; click a header to re-sort. ",
      "n = bins where both the feature and the rate are present. Sliding-window bins ",
      "overlap and are spatially autocorrelated, so p-values are not reported."),
      format(n, big.mark = ","), corr_scope_label()))
  })
  
  # Correlation matrix, shown in a pop-up window
  observeEvent(input$corr_open, {
    showModal(modalDialog(
      title = "Correlation matrix",
      size = "xl", easyClose = TRUE, footer = modalButton("Close"),
      fluidRow(
        column(3, radioButtons("cm_method", "Method",
                               c("Spearman" = "spearman", "Pearson" = "pearson"),
                               selected = input$cm_method %||% "spearman", inline = TRUE)),
        column(5, radioButtons("cm_vars", "Variables",
                               c("Recombination rate + all features" = "all",
                                 "Variables selected for plotting" = "ys"),
                               selected = input$cm_vars %||% "all", inline = TRUE)),
        column(4, radioButtons("cm_order", "Order",
                               c("Feature category" = "category",
                                 "Hierarchical clustering" = "cluster"),
                               selected = input$cm_order %||% "category", inline = TRUE))
      ),
      uiOutput("cm_note"),
      plotlyOutput("corr_matrix", height = "auto")
    ))
  })
  
  corr_mat <- reactive({
    req(input$cm_method, input$cm_vars, input$cm_order, input$corr_target)
    d <- corr_data(); fc <- feature_cols(); target <- input$corr_target
    vars <- if (input$cm_vars == "ys") unique(c(target, intersect(input$ys, names(d))))
    else c(target, fc$features)
    nm <- numeric_matrix(d, vars)
    validate(need(ncol(nm$M) >= 2,
                  "At least two variables with non-constant values are needed."))
    R <- withProgress(message = "Computing correlation matrix...", value = 0.5,
                      cor_matrix(nm$M, input$cm_method))
    N <- crossprod(!is.na(nm$M))                       # pairwise-complete bins per cell
    v <- colnames(R)
    ord <- if (input$cm_order == "cluster") cluster_order(R) else {
      grp <- ifelse(v %in% fc$targets, 0L,
                    match(vapply(v, feature_category, character(1)), CATEGORY_LEVELS))
      order(grp)                                       # stable: target(s) first, then by category
    }
    list(R = R[ord, ord, drop = FALSE], N = N[ord, ord, drop = FALSE],
         target = target, dropped = nm$dropped, n = nrow(d))
  })
  
  output$cm_note <- renderUI({
    cm <- corr_mat()
    helpText(sprintf(paste0(
      "%s correlation over %s bins (%s), using the bins where both variables are present ",
      "(hover a cell for n). Bins and recombination rate are set in the correlation panel.%s"),
      c(spearman = "Spearman", pearson = "Pearson")[[input$cm_method]],
      format(cm$n, big.mark = ","), corr_scope_label(),
      if (length(cm$dropped))
        paste0(" Left out (constant or < 3 values): ", paste(cm$dropped, collapse = ", "), ".")
      else ""))
  })
  
  output$corr_matrix <- renderPlotly({
    cm <- corr_mat(); R <- cm$R; v <- colnames(R); p <- length(v)
    sym <- if (input$cm_method == "spearman") "\u03c1" else "r"
    hover <- matrix(sprintf("%s \u00d7 %s<br>%s = %s<br>n = %s",
                            rep(v, times = p), rep(v, each = p), sym,
                            ifelse(is.na(R), "NA", sprintf("%.3f", R)),
                            format(cm$N, big.mark = ",", trim = TRUE)), nrow = p)
    # Grow the plot with the number of variables so every label stays readable
    fs <- if (p > 60) 9 else if (p > 35) 10 else 12
    h  <- max(650, 11 * p + 220)
    ticks <- ifelse(v == cm$target, paste0("<b>", v, "</b>"), v)   # highlight the rate
    axis <- list(title = "", type = "category", categoryorder = "array", categoryarray = v,
                 tickvals = v, ticktext = ticks, tickfont = list(size = fs),
                 showgrid = FALSE, zeroline = FALSE, automargin = TRUE, constrain = "domain")
    # Print the coefficient in each cell when the matrix is small enough to read
    anns <- list()
    if (p <= 20) {
      idx <- which(!is.na(R), arr.ind = TRUE)
      anns <- lapply(seq_len(nrow(idx)), function(k) {
        i <- idx[k, 1]; j <- idx[k, 2]
        list(x = j - 1, y = i - 1, text = sprintf("%.2f", R[i, j]), showarrow = FALSE,
             font = list(size = 10, color = if (abs(R[i, j]) > 0.6) "#ffffff" else "#333333"))
      })
    }
    plot_ly(x = v, y = v, z = R, text = hover, type = "heatmap", hoverinfo = "text",
            zmin = -1, zmax = 1, xgap = 1, ygap = 1, height = h,
            colorscale = list(c(0, CORR_POLES[["neg"]]), c(0.5, CORR_POLES[["mid"]]),
                              c(1, CORR_POLES[["pos"]])),
            colorbar = list(title = sym, len = 0.6, tickvals = c(-1, -0.5, 0, 0.5, 1))) %>%
      layout(xaxis = c(axis, list(tickangle = if (p > 30) -90 else -45)),
             yaxis = c(axis, list(autorange = "reversed", scaleanchor = "x")),
             annotations = anns, plot_bgcolor = "#ffffff",
             margin = list(l = 10, r = 10, t = 10, b = 10)) %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(filename = "correlation_matrix", scale = 2))
  })
  

  # ---- Distributions & pairwise correlations of the plotted variables ----
  # Second tab: boxplots of each selected variable + scatterplot matrix.
  # The button in the left panel switches to the distributions tab
  observeEvent(input$dist_open, updateTabsetPanel(session, "main_tabs", selected = "dist"))

  # Bins for the selected donor / sex / dataset; current chromosome or all.
  # Recombinogenic bins = bins with >= 1 CO (Number_of_COs, or rate > 0 if the
  # CO count is missing).
  dv_data <- reactive({
    req(input$donor, input$male, input$dataset, input$dv_scope)
    d <- subset_rows(raw(), input$donor, input$male, input$dataset)
    if (input$dv_scope == "chr") {
      req(input$chr)
      d <- d[Chr_ID == input$chr]
    }
    vars <- intersect(input$ys, names(d))
    vars <- vars[vapply(vars, function(v) is.numeric(d[[v]]), logical(1))]
    validate(need(length(vars) > 0, "Select at least one variable to plot."),
             need(nrow(d) > 0, "No bins for this combination of selections."))
    rate <- intersect(vars, TARGET_COLS)
    rec <- if ("Number_of_COs" %in% names(d)) is.finite(d$Number_of_COs) & d$Number_of_COs > 0
           else if (length(rate)) is.finite(d[[rate[1]]]) & d[[rate[1]]] > 0
           else rep(FALSE, nrow(d))
    d <- copy(d)
    rec_vars <- paste0(rate, REC_SUFFIX)
    for (k in seq_along(rate))
      set(d, j = rec_vars[k], value = ifelse(rec, as.numeric(d[[rate[k]]]), NA_real_))
    # Boxplots: each rate followed by its recombinogenic-bin version
    box_vars <- vars
    for (k in seq_along(rate))
      box_vars <- append(box_vars, rec_vars[k], after = match(rate[k], box_vars))
    # Optional split of the bins by a feature threshold
    g <- NULL; labs <- NULL
    f <- input$dv_split
    if (isTRUE(nzchar(f) && f != "none" && f %in% names(d)) && is.numeric(input$dv_split_thr)) {
      thr <- input$dv_split_thr
      g <- ifelse(is.finite(d[[f]]), d[[f]] > thr, NA)
      labs <- c(sprintf("%s > %s", f, signif(thr, 4)), sprintf("%s \u2264 %s", f, signif(thr, 4)))
    }
    list(d = d, vars = vars, rec = rec, rate = rate, rec_vars = rec_vars,
         box_vars = box_vars, g = g, labs = labs)
  })

  dv_scope_label <- reactive({
    sex <- c("0" = "female meiosis", "1" = "male meiosis", "NA" = "both sexes")[[input$male]]
    where <- if (input$dv_scope == "chr") input$chr else "all chromosomes"
    paste(where, input$donor, sex, input$dataset, sep = " | ")
  })

  # Correlation matrix + data for the scatterplot matrix; rate(s) placed last so
  # the bottom row plots the rate (y) against every feature (x)
  dv_cor_for <- function(d, vars) {
    req(input$dv_method, input$dv_order)
    nm <- numeric_matrix(d, vars)
    validate(need(ncol(nm$M) >= 2,
                  "At least two plotted variables with non-constant values are needed for correlations."))
    R <- cor_matrix(nm$M, input$dv_method)
    ord <- if (input$dv_order == "cluster") cluster_order(R) else seq_len(ncol(R))
    is_t <- is_rate_col(colnames(R)[ord])
    ord <- c(ord[!is_t], ord[is_t])
    list(R = R[ord, ord, drop = FALSE], M = nm$M[, ord, drop = FALSE],
         dropped = nm$dropped, n = nrow(d))
  }

  # All bins, with the plain rate
  dv_cor <- reactive({
    x <- dv_data()
    c(dv_cor_for(x$d, x$vars), list(g = x$g, labs = x$labs))
  })

  # Recombinogenic bins only, with the rate renamed *_recombinogenic
  dv_cor_rec <- reactive({
    x <- dv_data()
    validate(need(sum(x$rec) >= 3, "Fewer than 3 recombinogenic bins (bins with ≥ 1 CO)."))
    vars <- x$vars
    vars[match(x$rate, vars)] <- x$rec_vars
    c(dv_cor_for(x$d[x$rec], vars),
      list(g = if (is.null(x$g)) NULL else x$g[x$rec], labs = x$labs))
  })

  output$dv_note <- renderUI({
    x <- dv_data()
    helpText(sprintf(paste0(
      "%s bins (%s), of which %s (%.1f%%) are recombinogenic (≥ 1 CO). Variables are ",
      "those selected for plotting. Boxplots: one panel per variable with its own y scale ",
      "(box = IQR, line = median, whiskers = 1.5 × IQR); features use all bins, and ",
      "%s shows the rate in recombinogenic bins only. Matrices: scatter plots with a ",
      "lowess curve (red) below the diagonal (up to %s bins drawn at random), densities ",
      "on the diagonal, and coefficients above it, computed on all bins of the matrix ",
      "where both variables are present. Stars: * p < 0.05, ** p < 0.01, *** p < 0.001. ",
      "Bins are spatially autocorrelated (and sliding windows overlap), so the stars ",
      "overstate significance; read them as a rough guide."),
      format(nrow(x$d), big.mark = ","), dv_scope_label(),
      format(sum(x$rec), big.mark = ","), 100 * mean(x$rec),
      if (length(x$rec_vars)) paste(x$rec_vars, collapse = " / ") else "the *_recombinogenic rate",
      format(PAIRS_MAX_POINTS, big.mark = ",")))
  })

  dropped_note <- function(cm, what) {
    txt <- sprintf("%s bins (%s).", format(cm$n, big.mark = ","), what)
    if (length(cm$dropped))
      txt <- paste0(txt, " Left out (constant or < 3 values): ",
                    paste(cm$dropped, collapse = ", "), ".")
    helpText(txt)
  }
  output$dv_pairs_note     <- renderUI(dropped_note(dv_cor(), "all bins"))
  output$dv_pairs_rec_note <- renderUI(dropped_note(dv_cor_rec(), "bins with ≥ 1 CO"))

  dv_box_ncol <- reactive(min(6, length(dv_data()$box_vars)))

  output$dv_box <- renderPlot({
    x <- dv_data(); bv <- x$box_vars
    w <- x$d[, bv, with = FALSE]
    split <- !is.null(x$g)
    if (split) w[, .grp := factor(ifelse(x$g, x$labs[1], x$labs[2]), levels = x$labs)]
    long <- melt(w, measure.vars = bv, id.vars = if (split) ".grp" else character(0),
                 variable.name = "variable", value.name = "value", variable.factor = FALSE)
    long <- long[is.finite(value)]
    if (split) long <- long[!is.na(.grp)]
    validate(need(nrow(long) > 0, "No finite values to plot."))
    long[, variable := factor(variable, levels = bv)]
    long[, kind := ifelse(variable %in% x$rec_vars, "Recombination rate (recombinogenic bins)",
                   ifelse(is_rate_col(variable), "Recombination rate (all bins)", "Feature"))]
    base <- if (split)
      ggplot(long, aes(x = .grp, y = value, fill = .grp)) +
        scale_fill_manual(values = setNames(unname(SPLIT_COLORS), x$labs), name = NULL)
    else
      ggplot(long, aes(x = "", y = value, fill = kind)) +
        scale_fill_manual(values = c("Recombination rate (all bins)" = "#e9c46a",
                                     "Recombination rate (recombinogenic bins)" = "#f4a261",
                                     "Feature" = "grey80"), name = NULL)
    # n per box, printed above it
    cnt <- long[, .(n = .N), by = c("variable", if (split) ".grp" else "kind")]
    cnt[, lab := paste0("n = ", format(n, big.mark = ",", trim = TRUE))]
    base +
      geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4, width = 0.6) +
      geom_text(data = cnt, vjust = 1.3, size = 3.4, colour = "grey25", inherit.aes = FALSE,
                mapping = if (split) aes(x = .grp, y = Inf, label = lab)
                          else aes(x = "", y = Inf, label = lab)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
      facet_wrap(~variable, scales = "free_y", ncol = dv_box_ncol()) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top", panel.grid.major.x = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_text(face = "bold"))
  }, height = function() 60 + 230 * ceiling(length(dv_data()$box_vars) / dv_box_ncol()))

  # Square scatterplot matrices: as tall as the tab is wide
  pairs_height <- function(id, cm_fun) function() {
    w <- session$clientData[[paste0("output_", id, "_width")]]
    p <- tryCatch(ncol(cm_fun()$R), error = function(e) 2)
    if (is.null(w) || w <= 0) max(500, 95 * p + 60) else max(500, w)
  }

  output$dv_pairs <- renderPlot({
    cm <- dv_cor()
    withProgress(message = "Drawing correlation matrix (all bins)...", value = 0.5,
                 draw_pairs(cm$M, cm$R, input$dv_method, cm$g, cm$labs))
  }, res = 110, height = pairs_height("dv_pairs", dv_cor))

  output$dv_pairs_rec <- renderPlot({
    cm <- dv_cor_rec()
    withProgress(message = "Drawing correlation matrix (recombinogenic bins)...", value = 0.5,
                 draw_pairs(cm$M, cm$R, input$dv_method, cm$g, cm$labs))
  }, res = 110, height = pairs_height("dv_pairs_rec", dv_cor_rec))

  # ---- Feature vs. CO state (Poisson test) ----
  # Slider range follows the selected feature; default threshold reset on change
  # Threshold slider: range follows the feature; value reset to the default
  update_thr_slider <- function(id, f) {
    req(f, f %in% names(raw()), is.numeric(raw()[[f]]))
    x <- raw()[[f]]; x <- x[is.finite(x)]
    req(length(x) > 0)
    r <- range(x); if (diff(r) == 0) r <- r + c(-0.5, 0.5)
    step <- 10^floor(log10(diff(r) / 100))
    lo <- floor(r[1] / step) * step; hi <- ceiling(r[2] / step) * step
    thr <- min(max(default_threshold(f, x), lo), hi)
    updateSliderInput(session, id, min = lo, max = hi, step = step,
                      value = round(thr / step) * step)
  }
  observeEvent(list(input$poi_feature, raw()), update_thr_slider("poi_thr", input$poi_feature))
  observeEvent(list(input$dv_split, raw()), {
    req(input$dv_split != "none")
    update_thr_slider("dv_split_thr", input$dv_split)
  })
  
  # Non-overlapping bins for the selected donor / sex / dataset, all chromosomes
  poi_bins <- reactive({
    req(input$donor, input$male, input$dataset, input$poi_state)
    d <- subset_rows(raw(), input$donor, input$male, input$dataset)
    validate(need("Number_of_COs" %in% names(d), "Number_of_COs column not found."),
             need(nrow(d) > 0, "No rows for this combination of selections."))
    setorder(d, Chr_ID, Interval_start)
    if (input$poi_state %in% c("cold", "hot"))   # spot calls use the sliding-window bins
      d[, Spot_class := call_spots(.SD, input$min_len, input$hot_p, input$cold_p)$class,
        by = Chr_ID]
    d[d[, .I[tile_keep(Interval_start, Interval_end)], by = Chr_ID]$V1]
  })
  
  poi_res <- reactive({
    f <- input$poi_feature; thr <- input$poi_thr
    req(f, is.numeric(thr), input$poi_scope)
    d <- poi_bins()
    validate(need(f %in% names(d), "Feature not found."))
    keep <- is.finite(d[[f]]) & is.finite(d$Number_of_COs) & !is.na(d$Chr_ID) &
      is.finite(d$Interval_start) & is.finite(d$Interval_end)
    d <- d[keep]
    validate(need(nrow(d) > 0, "No bins with values for this feature."))
    pos <- d[[f]] > thr
    state <- switch(input$poi_state,
                    noco = d$Number_of_COs == 0, co = d$Number_of_COs > 0,
                    cold = d$Spot_class %in% "Coldspot", hot = d$Spot_class %in% "Hotspot")
    state[is.na(state)] <- FALSE
    len <- d$Interval_end - d$Interval_start
    co  <- d$Number_of_COs
    chr <- as.character(d$Chr_ID)
    rt  <- intersect(c(input$corr_target, TARGET_COLS), names(d))[1]
    rate <- if (is.na(rt)) NULL else d[[rt]]
    if (input$poi_scope == "all")
      return(cbind(Chromosome = "All", feature_state_test(co, len, pos, state, chr, rate)))
    do.call(rbind, lapply(unique(chr), function(cc) {
      i <- chr == cc
      cbind(Chromosome = cc, feature_state_test(co[i], len[i], pos[i], state[i], chr[i],
                                                rate[i]))
    }))
  })
  
  # Bonferroni: each p-value column is multiplied by the number of rows tested
  # (chromosomes, or 1 when pooled) x the number of features screened, capped at 1
  POI_P_COLS <- c("Poisson_p", "Fisher_p", "QP_p", "Wilcox_p")
  POI_LABELS <- c(
    Chromosome = "Chr", Bins_pos = "Bins with feature", Bins = "Bins",
    O = "Observed COs", E = "Expected COs", OE = "O/E", Poisson_p = "Poisson p",
    Poisson_p_bonf = "Poisson p (Bonf.)", Precision = "Precision", Recall = "Recall",
    OR = "Odds ratio", Fisher_p = "Fisher p", Fisher_p_bonf = "Fisher p (Bonf.)",
    Rate_pos = "Mean rate (feature)", Rate_neg = "Mean rate (no feature)",
    QP_RR = "Rate ratio (quasi-Poisson)", QP_p = "Quasi-Poisson p",
    QP_p_bonf = "Quasi-Poisson p (Bonf.)",
    RecRate_pos = "Mean rate, ≥ 1 CO (feature)",
    RecRate_neg = "Mean rate, ≥ 1 CO (no feature)",
    Wilcox_p = "Wilcoxon p, ≥ 1 CO", Wilcox_p_bonf = "Wilcoxon p, ≥ 1 CO (Bonf.)")

  poi_n_tests <- reactive({
    k <- input$poi_ntests
    k <- if (is.numeric(k) && is.finite(k) && k >= 1) floor(k) else 1
    nrow(poi_res()) * k
  })

  output$poi_table <- renderDT({
    out <- as.data.frame(poi_res()); m <- poi_n_tests()
    for (pc in POI_P_COLS) {                       # adjusted column right after each p
      adj <- pmin(1, out[[pc]] * m)
      at <- match(pc, names(out))
      out <- cbind(out[seq_len(at)], setNames(data.frame(adj), paste0(pc, "_bonf")),
                   out[-seq_len(at)])
    }
    idx <- function(cols) match(intersect(cols, names(out)), names(out)) - 1   # 0-based
    datatable(out, rownames = FALSE, selection = "none", class = "compact hover",
              colnames = unname(POI_LABELS[names(out)]),
              options = list(dom = "t", paging = FALSE, scrollX = TRUE, scrollY = "260px",
                             scrollCollapse = TRUE, order = list(),
                             columnDefs = list(
                               list(targets = idx(c("Bins_pos", "Bins", "O")), render = dt_num(0)),
                               list(targets = idx("E"), render = dt_num(1)),
                               list(targets = idx(c("OE", "Precision", "Recall", "Rate_pos",
                                                    "Rate_neg", "QP_RR", "RecRate_pos",
                                                    "RecRate_neg")), render = dt_num(3)),
                               list(targets = idx(c("OR", POI_P_COLS, paste0(POI_P_COLS, "_bonf"))),
                                    render = dt_num(signif = TRUE))))) %>%
      formatStyle(paste0(POI_P_COLS, "_bonf"), fontWeight = "bold")
  })
  
  output$poi_note <- renderUI({
    req(poi_res())
    st <- c(noco = "no CO", co = "\u2265 1 CO", cold = "coldspot", hot = "hotspot")[[input$poi_state]]
    sex <- c("0" = "female meiosis", "1" = "male meiosis", "NA" = "both sexes")[[input$male]]
    helpText(sprintf(paste0(
      "%s | %s | %s. Non-overlapping bins; feature present if %s > %s. ",
      "Expected COs: each chromosome's COs spread evenly along its length; the Poisson test ",
      "compares observed with expected COs in bins with the feature (O/E < 1 = depleted). ",
      "Precision = P(%s | feature), recall = P(feature | %s); odds ratio with Fisher's exact ",
      "test (0.5 added to all cells if one is zero).%s ",
      "Rate columns use %s: mean in bins with / without the feature; rate ratio and p-value ",
      "from a quasi-Poisson model of COs (offset: bin length; per chromosome when pooled), ",
      "which allows for the extra variance from the many zero bins and hotspots and is ",
      "therefore more conservative than the Poisson test; and the same comparison within ",
      "recombinogenic bins (\u2265 1 CO), with a Wilcoxon rank-sum test. Neighbouring bins ",
      "are not independent, so treat all p-values as approximate. ",
      "Bonf. = Bonferroni-adjusted p (\u00d7 %s tests: %s row(s) \u00d7 %s feature(s) screened, capped at 1)."),
      input$donor, sex, input$dataset, input$poi_feature, signif(input$poi_thr, 4), st, st,
      if (input$poi_state %in% c("cold", "hot"))
        " Hot/coldspots use the hotspot / coldspot settings in the controls panel." else "",
      intersect(c(input$corr_target, TARGET_COLS), names(poi_bins()))[1],
      poi_n_tests(), nrow(poi_res()), poi_n_tests() / nrow(poi_res())))
  })
  
  # ---- Table ----
  output$sel_info <- renderText({
    d <- dat()$d; n <- length(selected())
    if (n) sprintf("Showing %d selected bin(s) out of %d.", n, nrow(d))
    else sprintf("No points selected: showing all %d bins for this selection.", nrow(d))
  })
  
  output$table <- renderDT({
    d <- dat()$d; sel <- selected()
    out <- if (length(sel)) d[.id %in% sel] else d
    out <- out[, setdiff(names(out), ".id"), with = FALSE]
    num <- names(out)[vapply(out, is.numeric, logical(1))]
    ints <- c("Interval_start", "Interval_end", "Interval_midpoint", "Bin_size_bp",
              "Number_of_COs", "Number_of_hybrids", "Male")
    dt <- datatable(out, rownames = FALSE, selection = "none",
                    options = list(scrollX = TRUE, pageLength = 15))
    fmt <- setdiff(num, ints)
    if (length(fmt)) dt <- formatSignif(dt, fmt, 5)
    dt
  })
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

shinyApp(ui, server)