#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Add FORMAT/AF to MuSE VCFs.
#          Called by the PAPIO process MUSE; moved from an inline heredoc.
# Usage  : muse_add_af.py <in.vcf> <out.vcf>
#------------------------------------------------------------------------------
# MuSE emits AD as (ref,alt) but no decimal AF. Compute AF = alt/(ref+alt).
import sys, pysam
vcf = pysam.VariantFile(sys.argv[1])
if "AF" not in vcf.header.formats:
    vcf.header.formats.add("AF", "1", "Float", "Allele frequency; computed from MuSE AD (alt/(ref+alt))")
out = pysam.VariantFile(sys.argv[2], "w", header=vcf.header)
for rec in vcf:
    for sn in rec.samples:
        s = rec.samples[sn]
        af = 0.0
        try:
            if "AD" in s and s["AD"] is not None:
                ad = s["AD"]
                if isinstance(ad, (tuple, list)) and len(ad) >= 2:
                    ref_d, alt_d = ad[0] or 0, ad[1] or 0
                    if (ref_d + alt_d) > 0:
                        af = round(alt_d / (ref_d + alt_d), 4)
        except Exception:
            pass
        s["AF"] = af
    out.write(rec)
out.close()
