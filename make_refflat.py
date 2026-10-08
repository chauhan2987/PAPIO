# GENCODE GTF -> Picard refFlat; exons written in ascending order (required by Picard, also for minus-strand genes)
import re, sys
from collections import defaultdict
gtf, out_path = sys.argv[1], sys.argv[2]
TID = re.compile(r'transcript_id "([^"]+)"'); GNM = re.compile(r'gene_name "([^"]+)"')
gene, chrom, strand, exons, cds = {}, {}, {}, defaultdict(list), {}
for line in open(gtf):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('exon', 'CDS'): continue
    tid = TID.search(f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
    if f[2] == 'exon':
        g = GNM.search(f[8]); gene[tid] = g.group(1) if g else tid
        chrom[tid], strand[tid] = f[0], f[6]; exons[tid].append((s, e))
    else:
        c0, c1 = cds.get(tid, (s, e)); cds[tid] = (min(c0, s), max(c1, e))
with open(out_path, 'w') as out:
    for tid in gene:
        ex = sorted(exons[tid]); tx0, tx1 = ex[0][0], ex[-1][1]
        c0, c1 = cds.get(tid, (tx1, tx1))                  # non-coding: CDS start = end = tx end
        out.write('\t'.join([gene[tid], tid, chrom[tid], strand[tid], str(tx0), str(tx1), str(c0), str(c1),
                             str(len(ex)), ','.join(str(a) for a, _ in ex) + ',', ','.join(str(b) for _, b in ex) + ',']) + '\n')
