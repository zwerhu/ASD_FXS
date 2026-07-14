#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
take  ES/SE/LP/EZ from the 10th column of GWAS VCF (FORMAT)to be scalars, and output as gzipped TSV。
- pvalue: p = 10^(-LP)
- zscore: z = ES/SE(if ES/SE then use EZ)
- effect_allele = ALT, non_effect_allele = REF
"""

import sys, gzip
from cyvcf2 import VCF

def get_scalar(rec, key):
    try:
        arr = rec.format(key)
        if arr is None:
            return None
        val = arr[0]
        if hasattr(val, "__len__"):
            val = val[0]
        return None if val is None else float(val)
    except Exception:
        return None

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.vcf.gz output.table.gz", file=sys.stderr)
        sys.exit(1)

    invcf, outgz = sys.argv[1], sys.argv[2]
    vcf = VCF(invcf)

    header = [
        "variant_id","chromosome","position",
        "effect_allele","non_effect_allele",
        "pvalue","zscore","effect_size","standard_error"
    ]
    with gzip.open(outgz, "wt") as out:
        out.write("\t".join(header) + "\n")
        for rec in vcf:
            chrom = str(rec.CHROM)
            pos = int(rec.POS)
            ref = rec.REF
            alt = rec.ALT[0] if rec.ALT else "NA"

            # variant_id： VCF ID；or chr:pos:ref:alt
            vid = rec.ID if rec.ID else f"chr{chrom}:{pos}:{ref}:{alt}"

            # FORMAT 
            es = get_scalar(rec, "ES")
            se = get_scalar(rec, "SE")
            lp = get_scalar(rec, "LP")
            ez = get_scalar(rec, "EZ")

            # calculate p、z
            p = None
            if lp is not None:
                p = 10 ** (-lp)           # p = 10^(-LP)
            else:
                # try INFO--PVAL
                p = rec.INFO.get("PVAL", None)
                if p is not None:
                    try: p = float(p)
                    except: p = None

            z = None
            if es is not None and se is not None and se != 0:
                z = es / se               # z = ES/SE
            elif ez is not None:
                z = ez

            def fmt(x):
                return "NA" if x is None else (str(x) if isinstance(x, str) else f"{x:.6g}")

            row = [
                vid, chrom, str(pos),
                alt, ref,
                fmt(p), fmt(z), fmt(es), fmt(se)
            ]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()