#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Fit COSMIC SBS signatures to one sample's SNVs with SigProfilerAssignment.
#          Called by mutational_signatures.sh (PAPIO process MUTATIONAL_SIGNATURES).
# Usage  : fit_signatures.py <vcf_dir> <out_dir> <genome_build> <cosmic_version> <exome: True|False>
#------------------------------------------------------------------------------
import sys
from SigProfilerAssignment import Analyzer as Analyze

###### GLOBAL VARIABLES
VCF_DIR, OUT_DIR, GENOME, COSMIC, EXOME = sys.argv[1:6]

##### Main Script
Analyze.cosmic_fit(
    samples=VCF_DIR,
    output=OUT_DIR,
    input_type="vcf",
    context_type="96",
    genome_build=GENOME,
    cosmic_version=float(COSMIC),
    exome=(EXOME == "True"),
    make_plots=True,
    verbose=False,
)
