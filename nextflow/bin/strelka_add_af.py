#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Add FORMAT/AF to Strelka2 somatic VCFs.
#          Called by the PAPIO process STRELKA_TAG_AF; moved from an inline heredoc.
# Usage  : strelka_add_af.py <in.vcf> <out.vcf> <snv|indel>
#------------------------------------------------------------------------------
import sys
import pysam

def main():
    vcf_in, vcf_out, variant_type = sys.argv[1], sys.argv[2], sys.argv[3]
    assert variant_type in ("snv", "indel"), "variant_type must be 'snv' or 'indel'"

    vcf = pysam.VariantFile(vcf_in)
    vcf.header.formats.add("AD", "R", "Integer",
        "Allelic depths (ref,alt); tier1 counts derived from Strelka U/TAR-TIR fields")
    vcf.header.formats.add("AF", "A", "Float",
        "Allele frequency per ALT; tier1 alt/(tier1 ref+tier1 alts), derived per Strelka user guide")

    out = pysam.VariantFile(vcf_out, "w", header=vcf.header)
    base_field = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}

    n_ok, n_fallback = 0, 0
    for rec in vcf:
        ref = rec.ref.upper()
        alts = [a.upper() for a in rec.alts] if rec.alts else []
        if not alts:                       # no ALT allele: nothing to tag
            out.write(rec); continue

        for sample_name in rec.samples:
            s = rec.samples[sample_name]
            ref_t1, alt_t1, alt_all = 0, 0, []
            try:
                if variant_type == "snv":
                    ref_field = base_field.get(ref)
                    if ref_field and ref_field in s and s[ref_field] is not None:
                        ref_t1 = s[ref_field][0]
                    for a in alts:
                        f = base_field.get(a)
                        alt_all.append(s[f][0] if f and f in s and s[f] is not None else 0)
                    alt_t1 = alt_all[0] if alt_all else 0
                else:
                    if "TAR" in s and s["TAR"] is not None:
                        ref_t1 = s["TAR"][0]
                    if "TIR" in s and s["TIR"] is not None:
                        alt_t1 = s["TIR"][0]
                n_ok += 1
            except Exception:
                n_fallback += 1

            alt_counts = [int(alt_t1)] + [int(x) for x in alt_all[1:]]
            alt_counts += [0] * (len(alts) - len(alt_counts))      # one count per ALT, always
            denom = ref_t1 + sum(alt_counts)
            # Number=R: one depth per allele, so bcftools norm -m- splits AD correctly.
            s["AD"] = tuple([int(ref_t1)] + alt_counts)
            s["AF"] = tuple(round(a / denom, 4) if denom > 0 else 0.0 for a in alt_counts)

        out.write(rec)

    out.close()
    sys.stderr.write(f"strelka_add_af: {n_ok} records tagged OK, {n_fallback} fallback-to-zero\n")

if __name__ == "__main__":
    main()
