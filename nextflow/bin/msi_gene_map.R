#!/usr/bin/env Rscript
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Map MSIsensor-pro unstable sites to genes.
#          Called by the PAPIO process MSI_GENE_MAP; moved from an inline heredoc.
# Usage  : msi_gene_map.R <somatic_sites> <gene_bed> <out.tsv>
#------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# Text columns read as character: type guessing turns an all-"T" repeat-unit
# column into logical TRUE (mononucleotide T repeats are common).
sites <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE,
                    colClasses=c(chromosome="character", repeat_unit_bases="character"))
req <- c("chromosome","location","repeat_times","repeat_unit_bases")
if (!all(req %in% colnames(sites)))
    stop(sprintf("unstable-site table missing columns: %s", paste(setdiff(req, colnames(sites)), collapse=",")), call.=FALSE)
bed <- read.delim(args[2], sep="\t", header=FALSE, stringsAsFactors=FALSE)
if (ncol(bed) < 5) stop("gene BED needs 5 columns: chr,start,end,symbol,ENSG", call.=FALSE)
bed <- bed[, 1:5]; colnames(bed) <- c("chr","start","end","gene","gene_id")
norm_chr <- function(x) { x <- sub("^chr", "", as.character(x)); x[x %in% c("M","MT")] <- "MT"; x }
sites$chrom_out <- sites$chromosome                      # reported as in the reference
sites$chromosome <- norm_chr(sites$chromosome); bed$chr <- norm_chr(bed$chr)   # matching only
bed$start      <- bed$start + 1L
sites$location <- suppressWarnings(as.integer(sites$location))
have_diff <- "difference" %in% colnames(sites); have_fdr <- "FDR" %in% colnames(sites)
bbc <- split(seq_len(nrow(bed)), bed$chr)
# Hits collected per site in a list and bound once (growing vectors with c() is quadratic
# in MSI-H WGS samples with tens of thousands of unstable sites).
hits <- vector("list", nrow(sites))
for (i in seq_len(nrow(sites))) {
    idx <- bbc[[ sites$chromosome[i] ]]; if (is.null(idx) || is.na(sites$location[i])) next
    j <- idx[bed$start[idx] <= sites$location[i] & bed$end[idx] >= sites$location[i]]
    if (length(j)) hits[[i]] <- cbind(site = i, gene = j)
}
hm   <- do.call(rbind, c(list(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("site", "gene")))), hits))
si   <- hm[, "site"]; gsym <- bed$gene[hm[, "gene"]]; gid <- bed$gene_id[hm[, "gene"]]
# Every column is sized to length(si): with no unstable site in a gene (the
# usual MSS case) si is empty and a scalar column (NA / "SOMATIC_UNSTABLE")
# would make data.frame() fail on 0 vs 1 rows. Zero hits -> header-only TSV.
n <- length(si)
out <- data.frame(
    Gene_ID=gid, Hugo_Symbol=gsym, Chromosome=sites$chrom_out[si], Position=sites$location[si],
    Repeat_Unit=sites$repeat_unit_bases[si], Ref_Repeat_Times=sites$repeat_times[si],
    Repeat_Diff = if (have_diff) sites$difference[si] else rep(NA, n),
    FDR         = if (have_fdr)  sites$FDR[si]        else rep(NA, n),
    Status      = rep("SOMATIC_UNSTABLE", n), stringsAsFactors=FALSE)
write.table(out, args[3], sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Mapped %d unstable site(s) to %d gene(s) (%d site-gene rows)\n", length(unique(si)), length(unique(gid)), nrow(out)))
