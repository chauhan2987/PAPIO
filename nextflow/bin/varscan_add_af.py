#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Add FORMAT/AF to VarScan2 VCFs.
#          Called by the PAPIO process VARSCAN; moved from an inline heredoc.
# Usage  : varscan_add_af.py <in.vcf> <out.vcf>
#------------------------------------------------------------------------------
# VarScan: RD=ref depth, AD=alt depth (single value), FREQ="NN.NN%".
# Header AD is redeclared to Number=R by sed upstream so pysam accepts a
# (ref,alt) pair. We only WRITE values here; every read defaults None -> 0
# so no None ever reaches int().
import sys, pysam

def first_num(x):
    # unwrap tuple/list, coerce None -> 0
    if x is None:
        return 0
    if isinstance(x, (tuple, list)):
        x = next((v for v in x if v is not None), 0)
    return x if x is not None else 0

def main():
    vcf = pysam.VariantFile(sys.argv[1])
    if "AF" not in vcf.header.formats:
        vcf.header.formats.add("AF", "1", "Float",
            "Allele frequency (0-1); from VarScan FREQ/AD")
    out = pysam.VariantFile(sys.argv[2], "w", header=vcf.header)
    for rec in vcf:
        for sn in rec.samples:
            s = rec.samples[sn]
            ref_d, alt_d, af = 0, 0, 0.0
            try:
                ref_d = first_num(s.get("RD"))
                # VarScan AD = ALT count; after Number=R reinterpretation it may
                # be a tuple with None padding - take last non-None value.
                ad_raw = s.get("AD")
                if isinstance(ad_raw, (tuple, list)):
                    alt_d = next((v for v in reversed(ad_raw) if v is not None), 0)
                else:
                    alt_d = ad_raw if ad_raw is not None else 0

                fq = s.get("FREQ")
                if fq is not None:
                    fq = fq[0] if isinstance(fq, (tuple, list)) else fq
                    try:
                        af = round(float(str(fq).rstrip("%")) / 100.0, 4)
                    except (ValueError, TypeError):
                        af = 0.0
                elif (ref_d + alt_d) > 0:
                    af = round(alt_d / (ref_d + alt_d), 4)
            except Exception:
                pass
            s["AD"] = (int(ref_d), int(alt_d))
            s["AF"] = af
        out.write(rec)
    out.close()

if __name__ == "__main__":
    main()
