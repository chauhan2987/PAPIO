#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Trim vcf2maf output to the PAPIO MAF columns and add VAFs.
#          Called by the PAPIO process VCF2MAF; moved from an inline heredoc.
# Usage  : trim_maf.py <out.maf> <germline|somatic>   (reads full.maf)
#------------------------------------------------------------------------------
import csv, sys
if len(sys.argv) != 3 or sys.argv[2] not in ("germline", "somatic"):
    sys.exit("usage: trim_maf.py <out.maf> <germline|somatic>")
germline = sys.argv[2] == "germline"
if germline:
    want = ["Hugo_Symbol","Entrez_Gene_Id","Gene","dbSNP_RS","NCBI_Build","Chromosome",
            "Start_Position","End_Position","Strand","Variant_Classification","Variant_Type",
            "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
            "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","HGVSc","HGVSp","HGVSp_Short",
            "Transcript_ID","Exon_Number","t_GT","n_GT","n_depth","n_ref_count","n_alt_count",
            "Allele","Protein_position","Amino_acids","Codons","EXON","INTRON","Normal_VAF",
            "all_effects","Feature","FILTER","CLIN_SIG","ClinVar_CLNSIG","ClinVar_CLNREVSTAT","MAX_AF",
            "gnomAD_fafmax_faf95_max_joint"]
    computed = {"Normal_VAF"}
else:
    want = ["Hugo_Symbol","Gene","dbSNP_RS","NCBI_Build","Chromosome",
            "Start_Position","End_Position","Strand","Variant_Classification",
            "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
            "HGVSc","HGVSp","HGVSp_Short","Transcript_ID","Exon_Number",
            "t_GT","n_GT","t_depth","t_ref_count","t_alt_count",
            "n_depth","n_ref_count","n_alt_count","Allele","Protein_position",
            "Amino_acids","Codons","EXON","INTRON","Tumor_VAF","Normal_VAF",
            "all_effects","Feature","FILTER"]
    computed = {"Tumor_VAF", "Normal_VAF"}

with open("full.maf") as f:
    lines = f.readlines()
comments = [l for l in lines if l.startswith("#")]
data     = [l for l in lines if not l.startswith("#")]
reader = csv.DictReader(data, delimiter="\t")
have = set(reader.fieldnames or [])
miss = [c for c in want if c not in computed and c not in have]
if miss:
    sys.stderr.write("WARN missing columns: " + ",".join(miss) + "\n")
def vaf(a, r):
    # alt / (ref + alt) from the allele counts, as in FILTER_SOMATIC_VCF and the classification
    # table (vcf2maf's *_depth is FORMAT/DP, which also counts uninformative reads)
    try:
        a, r = float(a), float(r)
        return round(a / (a + r), 4) if a + r > 0 else ""
    except (ValueError, TypeError):
        return ""
with open(sys.argv[1], "w", newline="") as out:
    for c in comments:
        out.write(c)
    w = csv.DictWriter(out, fieldnames=want, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()

    for row in reader:
        r = {c: row.get(c, "") for c in want if c not in computed}
        if not germline:
            r["Tumor_VAF"] = vaf(row.get("t_alt_count"), row.get("t_ref_count"))
        r["Normal_VAF"] = vaf(row.get("n_alt_count"), row.get("n_ref_count"))
        w.writerow(r)
