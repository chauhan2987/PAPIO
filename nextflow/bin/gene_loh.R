#!/usr/bin/env Rscript
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Gene- and arm-level copy number / LOH from Sequenza segments.
#          Called by the PAPIO process SEQUENZA gene/arm summary; moved from an inline heredoc.
# Usage  : gene_loh.R <segments> <gene_bed> <out_gene.tsv> <sample> <out_summary.tsv> <arm_bed> <arm_min> <out_arm.tsv> <confints> <sex>
#------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
seg <- read.delim(args[1], sep="\t", header=TRUE, stringsAsFactors=FALSE)
req <- c("chromosome","start.pos","end.pos","CNt","A","B")
if (!all(req %in% colnames(seg))) stop(sprintf("segments missing: %s", paste(setdiff(req,colnames(seg)),collapse=",")), call.=FALSE)
bed <- read.delim(args[2], sep="\t", header=FALSE, stringsAsFactors=FALSE)
if (ncol(bed) < 5) stop("gene BED needs 5 cols: chr,start,end,symbol,ENSG", call.=FALSE)
bed <- bed[,1:5]; colnames(bed) <- c("chr","start","end","gene","gene_id")
nc <- function(x){ x <- sub("^chr","",as.character(x)); x[x %in% c("M","MT")] <- "MT"; x }
bed$chr_out <- as.character(bed$chr)                          # reported as in the reference
seg$chromosome <- nc(seg$chromosome); bed$chr <- nc(bed$chr)   # matching only
bed$start <- bed$start + 1L
for (cc in c("CNt","A","B","start.pos","end.pos")) seg[[cc]] <- suppressWarnings(as.numeric(seg[[cc]]))

# Tumour ploidy = length-weighted mean autosomal CNt; expected CN baseline = round(ploidy).
auto_re <- "^([1-9]|1[0-9]|2[0-2])$"
seglen  <- seg$end.pos - seg$start.pos + 1
ok      <- grepl(auto_re, seg$chromosome) & !is.na(seg$CNt)
ploidy  <- if (any(ok)) sum(seglen[ok] * seg$CNt[ok]) / sum(seglen[ok]) else NA_real_
base    <- if (is.na(ploidy)) 2 else floor(ploidy + 0.5)

clf <- function(cn,b){ if(is.na(cn)&&is.na(b))return("NO_COVERAGE"); if(is.na(b))return("UNKNOWN"); if(!is.na(cn)&&cn==0)return("HOMDEL"); if(b==0)return("LOH"); "HET" }
sbc <- split(seq_len(nrow(seg)), seg$chromosome)
rc <- rep(NA_real_, nrow(bed)); ra <- rc; rb <- rc
for (i in seq_len(nrow(bed))) {
  idx <- sbc[[ bed$chr[i] ]]; if (is.null(idx)) next
  ov <- pmin(bed$end[i], seg$end.pos[idx]) - pmax(bed$start[i], seg$start.pos[idx]) + 1
  hit <- ov > 0; if (!any(hit)) next
  j <- idx[hit][ which.max(ov[hit]) ]; rc[i] <- seg$CNt[j]; ra[i] <- seg$A[j]; rb[i] <- seg$B[j]
}
out <- data.frame(Gene_ID=bed$gene_id, Hugo_Symbol=bed$gene, Chromosome=bed$chr_out, Start=bed$start, End=bed$end,
                  Total_CN=rc, Major_Allele=ra, Minor_Allele=rb, stringsAsFactors=FALSE)
out$LOH_Status <- mapply(clf, out$Total_CN, out$Minor_Allele)
# CN-neutral relative to tumour ploidy (not CN==2): a tetraploid A=4,B=0 is CN-neutral LOH.
out$LOH_Detail <- ifelse(out$LOH_Status != "LOH", out$LOH_Status,
                   ifelse(out$Total_CN == base, "CN_neutral_LOH",
                   ifelse(out$Total_CN <  base, "LOH_with_loss", "LOH_with_gain")))
# Males: chrX/chrY are single-copy by nature - B==0 there is expected, not LOH.
# Unknown sex: X/Y copy number cannot be interpreted.
xy_loh <- bed$chr %in% c("X", "Y") & out$LOH_Status == "LOH"
if (identical(args[10], "1")) {
  out$LOH_Status[xy_loh] <- "HEMIZYGOUS"
  out$LOH_Detail[xy_loh] <- "HEMIZYGOUS_EXPECTED"
} else if (!identical(args[10], "0")) {
  out$LOH_Status[xy_loh] <- "SEX_UNKNOWN"
  out$LOH_Detail[xy_loh] <- "SEX_UNKNOWN"
}
write.table(out, args[3], sep="\t", quote=FALSE, row.names=FALSE)

# Genome-level summary: ploidy, fraction genome altered (CNt != baseline), LOH fraction (autosomes).
L    <- sum(seglen[ok])
fga  <- if (L > 0) sum(seglen[ok & seg$CNt != base]) / L else NA_real_
okb  <- ok & !is.na(seg$B)
floh <- if (sum(seglen[okb]) > 0) sum(seglen[okb & seg$B == 0]) / sum(seglen[okb]) else NA_real_
# Sequenza's own fit (confints_CP.txt, Sequenza 3.x: row 2 = point estimate).
cp <- tryCatch(read.delim(args[9], sep="\t", header=TRUE), error=function(e) NULL)
s_cell <- if (!is.null(cp) && nrow(cp) >= 2 && "cellularity" %in% colnames(cp)) cp$cellularity[2] else NA
s_ploi <- if (!is.null(cp) && nrow(cp) >= 2 && "ploidy.estimate" %in% colnames(cp)) cp$ploidy.estimate[2] else NA
summ <- data.frame(metric = c("sample_id","tumour_ploidy_weighted","cn_baseline","autosomal_bp_assessed",
                              "fraction_genome_altered","fraction_genome_loh",
                              "sequenza_cellularity","sequenza_ploidy_estimate"),
                   value  = c(args[4], round(ploidy, 3), base, format(L, scientific=FALSE), round(fga, 4), round(floh, 4),
                              s_cell, s_ploi))
write.table(summ, args[5], sep="\t", quote=FALSE, row.names=FALSE)

# Arm-level calls: fraction of each autosomal arm lost / gained vs baseline, and LOH fraction.
if (length(args) >= 6 && nzchar(args[6])) {
  minf <- as.numeric(args[7])
  arm  <- read.delim(args[6], sep="\t", header=FALSE, stringsAsFactors=FALSE, comment.char="#")
  if (ncol(arm) < 4) stop("arm BED needs 4 cols: chr,start,end,arm", call.=FALSE)
  arm <- arm[,1:4]; colnames(arm) <- c("chr","start","end","arm")
  arm$chr_out <- as.character(arm$chr); arm$chr <- nc(arm$chr); arm$start <- arm$start + 1L
  rows <- list()
  for (i in seq_len(nrow(arm))) {
    if (!grepl(auto_re, arm$chr[i])) next
    idx <- sbc[[ arm$chr[i] ]]; if (is.null(idx)) next
    ov <- pmin(arm$end[i], seg$end.pos[idx]) - pmax(arm$start[i], seg$start.pos[idx]) + 1
    h  <- ov > 0 & !is.na(seg$CNt[idx]); if (!any(h)) next
    o  <- ov[h]; cn <- seg$CNt[idx][h]; b <- seg$B[idx][h]; cov <- sum(o)
    lf <- sum(o[cn < base]) / cov; gf <- sum(o[cn > base]) / cov; hf <- sum(o[!is.na(b) & b == 0]) / cov
    call <- if (lf >= minf) "LOSS" else if (gf >= minf) "GAIN" else "NEUTRAL"
    rows[[length(rows) + 1]] <- data.frame(arm=arm$arm[i], chromosome=arm$chr_out[i], covered_bp=format(cov, scientific=FALSE),
                                           loss_frac=round(lf,4), gain_frac=round(gf,4), loh_frac=round(hf,4), call=call)
  }
  tab <- if (length(rows)) do.call(rbind, rows) else
         data.frame(arm=character(), chromosome=character(), covered_bp=numeric(), loss_frac=numeric(),
                    gain_frac=numeric(), loh_frac=numeric(), call=character())
  write.table(tab, args[8], sep="\t", quote=FALSE, row.names=FALSE)
}
