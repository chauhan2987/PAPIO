#!/usr/bin/env Rscript
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Sequenza fit (purity/ploidy/segments) from a binned seqz file.
#          Called by the PAPIO process SEQUENZA; moved from an inline heredoc.
# Usage  : run_sequenza.R <binned.seqz.gz> <sample_id> <assembly> <gamma> <kmin> <female TRUE|FALSE> [segment_filter_bp]
#------------------------------------------------------------------------------
suppressMessages(library(sequenza))
suppressMessages(library(iotools))

# sequenza 3.0.0 calls chunk.apply(parallel=) in four internal functions;
# iotools >=0.3 renamed that argument to CH.PARALLEL. Patch all four in
# sequenza's own namespace. Verified: complete set, all patch cleanly.
# No-op if iotools is old (still has 'parallel').
if ("CH.PARALLEL" %in% names(formals(iotools::chunk.apply))) {
    .targets <- c("gc.sample.stats", "read.seqz", "read.seqz.chr", "sequenza.extract")
    for (.t in .targets) {
        .orig <- get(.t, envir = asNamespace("sequenza"))
        .src  <- deparse(.orig)
        .src  <- gsub("parallel = parallel", "CH.PARALLEL = parallel", .src, fixed = TRUE)
        .fn   <- eval(parse(text = .src))
        environment(.fn) <- asNamespace("sequenza")
        assignInNamespace(.t, .fn, ns = "sequenza")
    }
    cat("Patched", length(.targets), "sequenza functions to use CH.PARALLEL.\n")
}

args      <- commandArgs(trailingOnly = TRUE)
seqz_file <- args[1]
prefix    <- args[2]
assembly  <- args[3]
gamma     <- as.integer(args[4])
kmin      <- as.integer(args[5])
if (length(args) < 6 || !(args[6] %in% c("TRUE", "FALSE")))
    stop("usage: run_sequenza.R <seqz> <sample_id> <assembly> <gamma> <kmin> <female TRUE|FALSE> [segment_filter_bp]")
female    <- as.logical(args[6])
# Minimum segment length used to fit purity/ploidy (Sequenza default 3e6 bp; lower it for
# small targeted panels, where few segments reach 3 Mb and the fit has no data).
seg_filter <- if (length(args) >= 7) as.numeric(args[7]) else 3e6

set.seed(1234)

seqz_data <- sequenza.extract(
    seqz_file,
    assembly = assembly,
    gamma    = gamma,
    kmin     = kmin,
    verbose  = TRUE
)

# sequenza.extract wraps copynumber segmentation in try(): an unsupported
# assembly silently yields ONE segment per chromosome. Stop instead.
nseg <- sapply(seqz_data$segments, nrow)
print(nseg)
if (all(nseg <= 1))
    stop("Sequenza produced one segment per chromosome: segmentation failed ",
         "(the copynumber package must support assembly '", assembly, "')")

# Sex-chromosome names as they appear in the data (chrX/chrY or X/Y): Sequenza's
# default XY = c("X","Y") does not match chrX/chrY, so chrY was fitted in female
# samples and a male chrX was treated as diploid.
xy_pfx <- if (any(grepl("^chr", seqz_data$chromosomes))) "chr" else ""
XY     <- c(X = paste0(xy_pfx, "X"), Y = paste0(xy_pfx, "Y"))
cat("Sex chromosomes:", XY, "\n")

CP <- sequenza.fit(
    seqz_data,
    cellularity = seq(0.1, 1, 0.01),
    ploidy      = seq(1, 7, 0.1),
    segment.filter = seg_filter,
    female      = female,
    XY          = XY
)

sequenza.results(
    sequenza.extract = seqz_data,
    cp.table         = CP,
    sample.id        = prefix,
    out.dir          = ".",
    female           = female,
    XY               = XY
)
