#!/usr/bin/env python3
"""
Create exon and intron BED6 files from Chen et al. 2020 ORF catalogue.

Dual-method approach:
  Method A (annotation): Ensembl GRCh38 REST API for ENST IDs;
                         GENCODE v49 overlap inference for TCONS IDs.
  Method B (SPALN):      Protein-to-genome spliced alignment against hg38.

Chromosome names in output use CM accessions (GRCh38) for pipeline compatibility.

Coordinate convention in Chen 2020 catalogue (0-based, BED-like):
  + strand: start < stop;  extraction hg19[chrom][start:stop]
  - strand: start > stop;  extraction revcomp(hg19[chrom][stop:start+1])
"""

import bisect
import os, sys, re, time, gzip, tempfile, subprocess, logging
import pandas as pd
import requests
from pyfaidx import Fasta

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
CATALOGUE    = "/pl/active/ShenLab_PL/haploid_genetics/Chen et al. 2020 catalogue.xlsx"
AA_SEQS      = "/pl/active/ShenLab_PL/haploid_genetics/Chen et al. 2020 AA seqs.xlsx"
CHAIN        = "/scratch/alpine/wist9668/reference/liftover/hg19ToHg38.over.chain"
HG19_FA      = "/scratch/alpine/wist9668/reference/hg19.fa"
HG38_FA      = "/pl/active/ShenLab_PL/haploid_genetics/reference/GCA_000001405.15_GRCh38_genomic.fna"
ASSEMBLY_RPT = "/pl/active/ShenLab_PL/haploid_genetics/reference/annotation/GCA_000001405.15_GRCh38_assembly_report.txt"
GENCODE_GTF  = "/scratch/alpine/wist9668/reference/gencode_v49/gencode.v49.annotation.gtf.gz"
OUT_DIR      = "/scratch/alpine/wist9668/reference/annotation/chen2020"
LIFTOVER_BIN = "/projects/wist9668/software/anaconda/envs/hap1/bin/liftOver"
SPALN_BIN    = "/projects/wist9668/software/anaconda/envs/hap1/bin/spaln"
FLANK        = 500
ENSEMBL_38   = "https://rest.ensembl.org"
ENSEMBL_37   = "https://grch37.rest.ensembl.org"

os.makedirs(OUT_DIR, exist_ok=True)

# ── Utilities ──────────────────────────────────────────────────────────────────
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq):
    return seq.translate(_COMP)[::-1]

def translate(seq):
    aas = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i+3], 'X')
        if aa == '*':
            break
        aas.append(aa)
    return ''.join(aas)

def build_chr_to_cm(report_path):
    """Map UCSC chr names (col 9) → CM accessions (col 4) from GRCh38 assembly report."""
    m = {}
    with open(report_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) > 9 and parts[9] != 'na':
                m[parts[9]] = parts[4]
    return m

# ── Step 1: Load and filter catalogue ─────────────────────────────────────────
log.info("Loading Chen et al. 2020 catalogue ...")
df_ipsc = pd.read_excel(CATALOGUE, sheet_name='iPSC', header=3)
df_hff  = pd.read_excel(CATALOGUE, sheet_name='HFF',  header=3)
df = pd.concat([df_ipsc, df_hff]).drop_duplicates(subset='ORF name').reset_index(drop=True)
df = df[(df['ORF type'] == 'new') & (df['amino acid length'] >= 10)].copy()
for col in ('start genomic coordinate', 'stop genomic coordinate', 'amino acid length'):
    df[col] = df[col].astype(int)

df['coding_nt']    = df['amino acid length'] * 3 + 3
df['genomic_span'] = (df['start genomic coordinate'] - df['stop genomic coordinate']).abs()
df['single_exon']  = df['genomic_span'] == df['coding_nt']
df['is_enst']      = df['transcript ID'].str.startswith('ENST')

log.info(f"  {len(df)} ORFs: {df['single_exon'].sum()} single-exon, "
         f"{(~df['single_exon']).sum()} spliced "
         f"(ENST={df['is_enst'].sum()}, TCONS={(~df['is_enst']).sum()})")

# ── Step 2: LiftOver hg19 → hg38 ──────────────────────────────────────────────
log.info("Running liftOver hg19 → hg38 ...")
hg38_spans = {}     # orf_name → {chrom_chr, start, end, strand}
unmapped   = set()

with tempfile.TemporaryDirectory() as tmp:
    hg19_bed    = os.path.join(tmp, 'hg19.bed')
    hg38_bed    = os.path.join(tmp, 'hg38.bed')
    unmap_bed   = os.path.join(tmp, 'unmapped.bed')

    with open(hg19_bed, 'w') as f:
        for _, row in df.iterrows():
            s = int(row['start genomic coordinate'])
            e = int(row['stop genomic coordinate'])
            strand = row['strand']
            # BED coords: 0-based half-open, chromStart < chromEnd
            if strand == '+':
                lo, hi = s, e          # start=0-based chromStart, stop=0-based chromEnd
            else:
                lo, hi = e, s + 1      # stop=0-based chromStart, start+1=0-based chromEnd
            f.write(f"{row['chromosome']}\t{lo}\t{hi}\t{row['ORF name']}\t0\t{strand}\n")

    ret = subprocess.run(
        [LIFTOVER_BIN, hg19_bed, CHAIN, hg38_bed, unmap_bed],
        capture_output=True, text=True
    )
    if ret.returncode not in (0, 1):
        log.error(f"liftOver failed: {ret.stderr}")
        sys.exit(1)

    with open(hg38_bed) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            hg38_spans[p[3]] = {
                'chrom_chr': p[0],
                'start': int(p[1]),   # 0-based
                'end':   int(p[2]),
                'strand': p[5],
            }

    with open(unmap_bed) as f:
        for line in f:
            if not line.startswith('#') and '\t' in line:
                unmapped.add(line.rstrip('\n').split('\t')[3])

log.info(f"  Mapped: {len(hg38_spans)}, unmapped: {len(unmapped)}")
df = df[df['ORF name'].isin(hg38_spans)].copy()

# ── Step 3: Chr→CM map ────────────────────────────────────────────────────────
chr_to_cm = build_chr_to_cm(ASSEMBLY_RPT)
log.info(f"  CM map loaded: {len(chr_to_cm)} entries")

# ── Step 3b: Load AA sequence catalogue ──────────────────────────────────────
log.info("Loading Chen et al. 2020 AA sequence catalogue ...")
aa_df = pd.read_excel(AA_SEQS)
orf_to_protein = dict(zip(aa_df['ORF name'], aa_df['Amino acid sequence']))
log.info(f"  {len(orf_to_protein)} ORF AA sequences loaded")

# ── Step 4: Load GENCODE v49 exons ────────────────────────────────────────────
log.info("Loading GENCODE v49 exons (this may take a minute) ...")
gencode_rows = []
with gzip.open(GENCODE_GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'exon':
            continue
        chrom  = p[0]
        start  = int(p[3]) - 1   # 0-based
        end    = int(p[4])
        strand = p[6]
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m:
            continue
        tx_id = m.group(1)
        gencode_rows.append((chrom, start, end, strand, tx_id))

gencode_df = pd.DataFrame(gencode_rows, columns=['chrom', 'start', 'end', 'strand', 'transcript_id'])
log.info(f"  {len(gencode_df)} exon rows, {gencode_df['transcript_id'].nunique()} transcripts")

# ── Step 5: Genome handles + Ensembl session ───────────────────────────────────
log.info("Opening genome FASTA handles ...")
hg19 = Fasta(HG19_FA, as_raw=True, sequence_always_upper=True)
hg38 = Fasta(HG38_FA, as_raw=True, sequence_always_upper=True)

session = requests.Session()
session.headers.update({'Content-Type': 'application/json'})

def ensembl_get(endpoint, server=ENSEMBL_38, retries=3):
    url = f"{server}/{endpoint}"
    for attempt in range(retries):
        try:
            r = session.get(url, headers={'Content-Type': 'application/json'}, timeout=30)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 429:
                time.sleep(float(r.headers.get('Retry-After', 5)))
            elif r.status_code == 400:
                return None
            else:
                time.sleep(2 ** attempt)
        except Exception as e:
            log.debug(f"API error ({url}): {e}")
            time.sleep(2 ** attempt)
    return None

def ensembl_seq(endpoint, server=ENSEMBL_38, retries=3):
    """Return plain-text sequence from Ensembl sequence endpoint."""
    url = f"{server}/{endpoint}"
    for attempt in range(retries):
        try:
            r = session.get(url, headers={'Content-Type': 'text/plain'}, timeout=30)
            if r.status_code == 200:
                return r.text.strip()
            elif r.status_code == 429:
                time.sleep(float(r.headers.get('Retry-After', 5)))
            elif r.status_code == 400:
                return None
            else:
                time.sleep(2 ** attempt)
        except Exception as e:
            log.debug(f"Seq API error ({url}): {e}")
            time.sleep(2 ** attempt)
    return None

# ── Method A helpers ──────────────────────────────────────────────────────────
def map_tx_to_genomic(exon_list_tx_order, strand, tx_start, tx_end):
    """
    Map ORF transcript positions [tx_start, tx_end) to genomic exon intervals.
    exon_list_tx_order: list of (g_start_0based, g_end_0based) in mRNA 5'→3' order.
    Returns list of (g_start, g_end) 0-based half-open intervals.
    """
    result = []
    tx_pos = 0
    for g_start, g_end in exon_list_tx_order:
        exon_len = g_end - g_start
        ov_s = max(tx_pos, tx_start)
        ov_e = min(tx_pos + exon_len, tx_end)
        if ov_s < ov_e:
            if strand == '+':
                result.append((g_start + (ov_s - tx_pos), g_start + (ov_e - tx_pos)))
            else:
                result.append((g_end - (ov_e - tx_pos), g_end - (ov_s - tx_pos)))
        tx_pos += exon_len
    return result

def method_a_enst(row, hg38_span):
    """Annotation exons for ENST ORFs via Ensembl GRCh38 REST API."""
    tx_id   = str(row['transcript ID'])
    tx_base = tx_id.split('.')[0]
    strand  = row['strand']
    tx_s    = int(row['start position on transcript'])
    tx_e    = int(row['stop position on transcript'])
    chrom_chr = hg38_span['chrom_chr']

    for server in [ENSEMBL_38, ENSEMBL_37]:
        data = ensembl_get(f"lookup/id/{tx_base}?expand=1&species=human", server=server)
        time.sleep(0.08)
        if data is None:
            continue

        exons = data.get('Exon', [])
        if not exons:
            continue

        seq_region = data.get('seq_region_name', '')
        if not seq_region.startswith('chr'):
            seq_region = 'chr' + seq_region

        # GRCh37 server → exons in hg19 coords; GRCh38 → hg38 coords
        # For GRCh37, we would need liftover of individual exons — complex.
        # Skip GRCh37 fallback and flag it; SPALN will still run.
        if server == ENSEMBL_37:
            return None, 'grch37_liftover_not_implemented'

        # Sort in transcript order (reverse for minus strand)
        raw = sorted([(int(e['start']) - 1, int(e['end'])) for e in exons],
                     key=lambda x: x[0] if strand == '+' else -x[0])

        orf_exons = map_tx_to_genomic(raw, strand, tx_s, tx_e)
        if not orf_exons:
            return None, 'tx_mapping_failed'

        return [(seq_region, gs, ge) for gs, ge in orf_exons], 'ensembl_grch38'

    return None, 'ensembl_api_failed'

def method_a_tcons(row, hg38_span):
    """Annotation exons for TCONS ORFs via GENCODE v49 overlap."""
    chrom_chr  = hg38_span['chrom_chr']
    span_s     = hg38_span['start']
    span_e     = hg38_span['end']
    strand     = row['strand']

    cands = gencode_df[
        (gencode_df['chrom']  == chrom_chr) &
        (gencode_df['strand'] == strand) &
        (gencode_df['start']  < span_e) &
        (gencode_df['end']    > span_s)
    ]
    if cands.empty:
        return None, 'no_gencode_overlap'

    # Score each transcript by total overlap with ORF span
    scores = {}
    for tx_id, grp in cands.groupby('transcript_id'):
        scores[tx_id] = sum(
            max(0, min(r['end'], span_e) - max(r['start'], span_s))
            for _, r in grp.iterrows()
        )

    best_tx = max(scores, key=scores.get)
    best_exons = gencode_df[gencode_df['transcript_id'] == best_tx].sort_values('start')

    result = []
    for _, exon in best_exons.iterrows():
        cs = max(exon['start'], span_s)
        ce = min(exon['end'],   span_e)
        if cs < ce:
            result.append((chrom_chr, cs, ce))

    if not result:
        return None, 'no_exons_after_clip'

    return result, f'GENCODE_v49:{best_tx}'

# ── Method B (SPALN) protein helpers ─────────────────────────────────────────
def protein_single_exon(row):
    """Exact protein for single-exon ORFs from hg19."""
    chrom  = row['chromosome']
    s      = int(row['start genomic coordinate'])
    e      = int(row['stop genomic coordinate'])
    strand = row['strand']
    try:
        if strand == '+':
            seq = str(hg19[chrom][s:e])
        else:
            seq = revcomp(str(hg19[chrom][e:s + 1]))
        return translate(seq), 'hg19_exact'
    except Exception as ex:
        log.debug(f"protein_single_exon failed {row['ORF name']}: {ex}")
        return None, 'hg19_extraction_failed'

def protein_enst_cds(tx_id_base):
    """Exact protein for ENST ORFs from Ensembl CDS sequence."""
    for server in [ENSEMBL_38, ENSEMBL_37]:
        cds = ensembl_seq(f"sequence/id/{tx_id_base}?type=cds", server=server)
        time.sleep(0.08)
        if cds:
            prot = translate(cds.upper())
            if prot:
                return prot, 'ensembl_cds'
    return None, 'ensembl_cds_failed'

def protein_tcons_bestfit(row):
    """
    Best-effort protein for spliced TCONS ORFs.
    Extract the genomic sequence from the annotated start position and translate
    in-frame, skipping stop codons (which likely mark intron boundaries).
    The resulting protein is imperfect near intron/exon junctions but sufficient
    to guide SPALN's GT-AG splice-site detection.
    """
    chrom  = row['chromosome']
    s      = int(row['start genomic coordinate'])
    e      = int(row['stop genomic coordinate'])
    strand = row['strand']
    aa_len = int(row['amino acid length'])
    # Extract a generous window around the ORF
    flank  = 300
    try:
        if strand == '+':
            lo = max(0, s - flank)
            seq = str(hg19[chrom][lo:e + flank])
            frame_start = s - lo   # offset of start codon in seq
        else:
            hi = s + 1
            lo = max(0, e - flank)
            seq = revcomp(str(hg19[chrom][lo:hi + flank]))
            frame_start = hi + flank - hi   # = flank = 0 offset after revcomp start

        # Translate in-frame from frame_start, skipping stop codons
        aas = []
        i = frame_start
        while i + 3 <= len(seq) and len(aas) < aa_len + 20:
            codon = seq[i:i + 3]
            aa = CODON_TABLE.get(codon, 'X')
            if aa != '*':
                aas.append(aa)
            i += 3

        if len(aas) >= max(5, aa_len // 2):
            return ''.join(aas[:aa_len]), 'hg19_bestfit'
    except Exception as ex:
        log.debug(f"protein_tcons_bestfit failed {row['ORF name']}: {ex}")
    return None, 'hg19_bestfit_failed'

def protein_from_aa_catalogue(orf_name):
    seq = orf_to_protein.get(orf_name)
    if seq:
        return str(seq), 'aa_catalogue'
    return None, 'not_in_catalogue'

# ── GT-AG search for TCONS multi-exon ORFs ───────────────────────────────────
def _gtag_validate(seq, exons, coding_nt, cat_protein):
    """Return True if concatenated exons form a valid ORF matching cat_protein."""
    cds = ''.join(seq[s:e] for s, e in exons)
    if len(cds) != coding_nt:
        return False
    if CODON_TABLE.get(cds[-3:], 'X') != '*':
        return False
    aas = [CODON_TABLE.get(cds[i:i+3], 'X') for i in range(0, coding_nt - 3, 3)]
    if '*' in aas:
        return False
    if cat_protein:
        return ''.join(aas) == cat_protein
    return True

def _check_partial(seq, exons_so_far, cat_protein):
    """Return False if completed codons in assembled exons contradict cat_protein."""
    if not cat_protein:
        return True
    cds = ''.join(seq[s:e] for s, e in exons_so_far)
    for i in range(len(cds) // 3):
        aa = CODON_TABLE.get(cds[i*3:(i+1)*3], 'X')
        if aa == '*' or aa != cat_protein[i]:
            return False
    return True

def _gtag_recurse(seq, gt_pos, ag_sorted, ag_set, cursor, introns_left, bp_left,
                  min_intron, coding_nt, cat_protein, exons_so_far, out, deadline):
    if len(out) > 10 or time.time() > deadline:
        return
    n = len(seq)
    if introns_left == 0:
        if bp_left != 0:
            return
        exons = exons_so_far + [(cursor, n)]
        if _gtag_validate(seq, exons, coding_nt, cat_protein):
            out.append(exons)
        return
    # Minimum space needed after this exon: remaining introns + 3 bp per remaining exon
    room_after = (introns_left - 1) * min_intron + introns_left * 3
    max_start  = n - bp_left - room_after
    for p in gt_pos:
        if p < cursor + 3:
            continue
        if p > max_start:
            break
        new_exons = exons_so_far + [(cursor, p)]
        if not _check_partial(seq, new_exons, cat_protein):
            continue
        if introns_left == 1:
            # Last intron: exact length bp_left from catalogue total_intron budget
            if bp_left < min_intron:
                continue
            q = p + bp_left
            if q > n - 3 or q not in ag_set:
                continue
            _gtag_recurse(seq, gt_pos, ag_sorted, ag_set, q, 0, 0,
                          min_intron, coding_nt, cat_protein, new_exons, out, deadline)
        else:
            # Variable intron: iterate only actual AG positions in valid budget window
            max_len = bp_left - (introns_left - 1) * min_intron
            lo = bisect.bisect_left(ag_sorted, p + min_intron)
            hi = bisect.bisect_right(ag_sorted, min(p + max_len, n - room_after))
            for q in ag_sorted[lo:hi]:
                _gtag_recurse(seq, gt_pos, ag_sorted, ag_set, q, introns_left - 1,
                              bp_left - (q - p), min_intron, coding_nt,
                              cat_protein, new_exons, out, deadline)
                if len(out) > 10 or time.time() > deadline:
                    return

def gtag_search_tcons(row, span_s, span_e, chrom_chr):
    """
    Find exact exon structure for a spliced TCONS ORF using GT-AG search in hg19.
    Uses the exact catalogue genomic span so non-ATG start codons work without SPALN.
    Only attempted when cat_protein is known (required for branch pruning).
    Returns list of (chrom_chr, hg38_start, hg38_end) or None.
    """
    chrom     = row['chromosome']
    s         = int(row['start genomic coordinate'])
    e         = int(row['stop genomic coordinate'])
    strand    = row['strand']
    orf       = row['ORF name']
    coding_nt = int(row['coding_nt'])

    cat_protein = orf_to_protein.get(orf)
    if cat_protein is None:
        return None  # no protein to validate against; too many false solutions

    try:
        if strand == '+':
            seq          = str(hg19[chrom][s:e]).upper()
            genomic_span = e - s
        else:
            seq          = revcomp(str(hg19[chrom][e:s])).upper()
            genomic_span = s - e
    except Exception:
        return None

    if len(seq) != genomic_span or genomic_span <= coding_nt:
        return None

    # total_intron is exact: derived from catalogue |start-stop| - coding_nt
    total_intron = genomic_span - coding_nt

    gt_pos    = [i for i in range(genomic_span - 1) if seq[i:i+2] == 'GT']
    ag_set    = {i + 2 for i in range(genomic_span - 1) if seq[i:i+2] == 'AG'}
    ag_sorted = sorted(ag_set)

    deadline = time.time() + 10  # safety net: abandon if still searching after 10s

    for min_intron in (60, 20):
        for n_introns in range(1, 5):
            solutions = []
            _gtag_recurse(seq, gt_pos, ag_sorted, ag_set, 0, n_introns, total_intron,
                          min_intron, coding_nt, cat_protein, [], solutions, deadline)
            if not solutions:
                continue

            if len(solutions) == 1:
                best = solutions[0]
            else:
                def _var(exons):
                    lens = [ex_e - ex_s for ex_s, ex_e in exons]
                    m = sum(lens) / len(lens)
                    return sum((l - m) ** 2 for l in lens)
                best = min(solutions, key=_var)

            # Map seq-space exon coordinates to hg38
            result = []
            if strand == '+':
                for ex_s, ex_e in best:
                    result.append((chrom_chr, span_s + ex_s, span_s + ex_e))
            else:
                # seq[k] corresponds to hg19 position s-1-k.
                # liftOver input was [e, s+1) → [span_s, span_e), length ≈ genomic_span+1.
                # Exon at seq[ex_s:ex_e] → hg38 [span_s+genomic_span-ex_e, span_s+genomic_span-ex_s)
                for ex_s, ex_e in best:
                    result.append((chrom_chr,
                                   span_s + genomic_span - ex_e,
                                   span_s + genomic_span - ex_s))
            log.debug(f"  gtag_search: {orf} solved with {n_introns} intron(s), "
                      f"min_intron={min_intron}, {len(solutions)} solution(s)")
            return result

    return None

# ── SPALN runner ──────────────────────────────────────────────────────────────
def run_spaln(protein, chrom_cm, span_s, span_e, strand, orf_name):
    """
    Align protein to hg38 genomic region via SPALN.
    Returns (list of (g_start_0based, g_end_0based), status_str) or (None, reason).
    """
    if not protein or len(protein) < 5:
        return None, 'protein_too_short'

    extract_s = max(0, span_s - FLANK)
    extract_e = span_e + FLANK

    try:
        genomic_seq = str(hg38[chrom_cm][extract_s:extract_e])
    except Exception as ex:
        return None, f'hg38_extract_failed:{ex}'

    if len(genomic_seq) < 30:
        return None, 'genomic_region_too_short'

    with tempfile.TemporaryDirectory() as tmp:
        prot_fa  = os.path.join(tmp, 'query.pep')
        g_fa     = os.path.join(tmp, 'target.fa')

        with open(prot_fa, 'w') as f:
            f.write(f'>{orf_name}\n{protein}\n')
        with open(g_fa, 'w') as f:
            f.write(f'>region_{orf_name}\n{genomic_seq}\n')

        spaln_s = '1' if strand == '+' else '2'
        # -Q0 uses DP mode; block-search modes (-Q4/-Q7) require larger genomes
        cmd = [SPALN_BIN, '-Q0', '-O0', f'-S{spaln_s}', g_fa, prot_fa]

        try:
            ret = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            return None, 'spaln_timeout'
        except Exception as ex:
            return None, f'spaln_error:{ex}'

        exons = parse_spaln_gff(ret.stdout, extract_s)

    if not exons:
        return None, 'spaln_no_exons'
    return exons, 'spaln_ok'

def parse_spaln_gff(gff_text, offset):
    """Parse SPALN -O0 GFF3 output → list of (chrom_0based_start, chrom_0based_end)."""
    exons = []
    for line in gff_text.splitlines():
        if not line or line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) < 9:
            continue
        if p[2].upper() not in ('CDS', 'EXON'):
            continue
        try:
            # GFF3 is 1-based inclusive → 0-based half-open
            gs = int(p[3]) - 1 + offset
            ge = int(p[4])     + offset
            exons.append((gs, ge))
        except ValueError:
            continue
    return sorted(set(exons))

# ── Comparison and confidence ─────────────────────────────────────────────────
def compare_exons(a, b, tol=3):
    """Compare two exon lists for boundary agreement (within tol bp)."""
    if a is None and b is None:
        return 'both_failed'
    if a is None:
        return 'annotation_failed'
    if b is None:
        return 'spaln_failed'
    if len(a) != len(b):
        return 'discordant_exon_count'
    for (as_, ae), (bs, be) in zip(
        sorted(a, key=lambda x: x[1]), sorted(b, key=lambda x: x[0])
    ):
        if abs(as_ - bs) > tol or abs(ae - be) > tol:
            return 'discordant'
    return 'high_confidence'

# ── Main processing loop ──────────────────────────────────────────────────────
log.info("Processing ORFs ...")
exon_bed_rows   = []
intron_bed_rows = []
report_rows     = []

total = len(df)
for idx, (_, row) in enumerate(df.iterrows()):
    orf = row['ORF name']
    if idx % 50 == 0:
        log.info(f"  {idx}/{total}: {orf}")

    span      = hg38_spans[orf]
    chrom_chr = span['chrom_chr']
    span_s    = span['start']
    span_e    = span['end']
    strand    = row['strand']
    is_single = bool(row['single_exon'])
    is_enst   = bool(row['is_enst'])
    chrom_cm  = chr_to_cm.get(chrom_chr)

    # ── Method A ────────────────────────────────────────────────────
    if is_single:
        method_a = [(chrom_chr, span_s, span_e)]
        annot_src = 'single_exon_liftover'
    elif is_enst:
        method_a, annot_src = method_a_enst(row, span)
    else:
        method_a, annot_src = method_a_tcons(row, span)

    # ── Method B (SPALN) ────────────────────────────────────────────
    if is_single:
        protein, prot_src = protein_single_exon(row)
    elif is_enst:
        tx_base = str(row['transcript ID']).split('.')[0]
        protein, prot_src = protein_enst_cds(tx_base)
        if protein is None:
            protein, prot_src = protein_from_aa_catalogue(orf)
        if protein is None:
            protein, prot_src = protein_tcons_bestfit(row)
            if prot_src == 'hg19_bestfit':
                prot_src = 'hg19_bestfit_enst_fallback'
    else:
        protein, prot_src = protein_from_aa_catalogue(orf)
        if protein is None:
            log.info(f"  No AA seq in catalogue for {orf}, using bestfit fallback")
            protein, prot_src = protein_tcons_bestfit(row)

    if chrom_cm is None:
        spaln_exons = None
        spaln_status = f'no_cm_for_{chrom_chr}'
    else:
        spaln_result, spaln_status = run_spaln(protein, chrom_cm, span_s, span_e, strand, orf)
        if spaln_result is not None:
            spaln_exons = [(chrom_chr, gs, ge) for gs, ge in spaln_result]
        else:
            spaln_exons = None

    # ── Compare, choose, output ─────────────────────────────────────
    # Strip chromosome from exon lists for comparison (coords only)
    a_coords = [(gs, ge) for (_, gs, ge) in method_a] if method_a else None
    b_coords = [(gs, ge) for (_, gs, ge) in spaln_exons] if spaln_exons else None

    confidence = compare_exons(a_coords, b_coords)

    # For TCONS multi-exon ORFs with exact AA seq: trust SPALN over GENCODE
    # (GENCODE exon structure does not match the original StringTie assembly).
    # For everything else: prefer Method A (annotation-anchored); fall back to SPALN.
    if not is_single and not is_enst and prot_src == 'aa_catalogue' and spaln_exons is not None:
        output_exons = spaln_exons
    elif confidence in ('high_confidence', 'discordant_exon_count', 'discordant'):
        output_exons = method_a     # use annotation; both available
    elif method_a is not None:
        output_exons = method_a
    elif spaln_exons is not None:
        output_exons = spaln_exons
    else:
        output_exons = None

    # ── CDS-length validation for ENST ORFs ──────────────────────────
    # If Ensembl exon structure gives wrong CDS length but SPALN (with exact
    # AA catalogue protein) gives the right length, prefer SPALN.
    if is_enst and not is_single and output_exons is not None and spaln_exons is not None:
        out_total   = sum(ge - gs for _, gs, ge in output_exons)
        spaln_total = sum(ge - gs for _, gs, ge in spaln_exons)
        expected    = int(row['coding_nt'])
        if abs(out_total - expected) > 3 and abs(spaln_total - expected) <= 3:
            output_exons = spaln_exons
            annot_src    = 'spaln_cds_validated'

    # ── GT-AG override for TCONS multi-exon ORFs ─────────────────────
    # Only invoked when SPALN CDS length doesn't match the catalogue coding_nt.
    # Uses exact hg19 coordinates + GT-AG splice site enumeration; handles
    # any start codon since we don't need SPALN to anchor to ATG.
    if not is_single and not is_enst:
        spaln_total = sum(ge - gs for _, gs, ge in spaln_exons) if spaln_exons else 0
        expected_cds = int(row['coding_nt'])
        if spaln_total != expected_cds:
            gtag = gtag_search_tcons(row, span_s, span_e, chrom_chr)
            if gtag is not None:
                output_exons = gtag
                annot_src    = 'gtag_hg19'
                confidence   = 'gtag_unique'

    # Write BED rows (CM chromosome)
    n_exons   = 0
    n_introns = 0
    if output_exons and chrom_cm:
        sorted_exons = sorted(output_exons, key=lambda x: x[1])
        for _, gs, ge in sorted_exons:
            exon_bed_rows.append(f"{chrom_cm}\t{gs}\t{ge}\t{orf}\t0\t{strand}")
            n_exons += 1
        for j in range(len(sorted_exons) - 1):
            intron_s = sorted_exons[j][2]
            intron_e = sorted_exons[j + 1][1]
            if intron_e > intron_s:
                intron_bed_rows.append(f"{chrom_cm}\t{intron_s}\t{intron_e}\t{orf}\t0\t{strand}")
                n_introns += 1

    report_rows.append({
        'orf_name':        orf,
        'single_exon':     is_single,
        'is_enst':         is_enst,
        'annotation_src':  annot_src,
        'protein_src':     prot_src,
        'spaln_status':    spaln_status,
        'confidence':      confidence,
        'exon_count':      n_exons,
        'intron_count':    n_introns,
    })

# ── Write outputs ─────────────────────────────────────────────────────────────
log.info("Writing output files ...")
exon_bed_path   = os.path.join(OUT_DIR, 'chen2020_exon_CM-mapped.bed')
intron_bed_path = os.path.join(OUT_DIR, 'chen2020_intron_CM-mapped.bed')
report_path     = os.path.join(OUT_DIR, 'chen2020_comparison_report.tsv')

with open(exon_bed_path, 'w') as f:
    f.write('\n'.join(exon_bed_rows) + '\n')
with open(intron_bed_path, 'w') as f:
    f.write('\n'.join(intron_bed_rows) + '\n')

report_df = pd.DataFrame(report_rows)
report_df.to_csv(report_path, sep='\t', index=False)

log.info(f"Exon BED:   {exon_bed_path}")
log.info(f"Intron BED: {intron_bed_path}")
log.info(f"Report:     {report_path}")
log.info("Summary:")
for k, v in report_df['confidence'].value_counts().items():
    log.info(f"  {k}: {v}")
log.info(f"  Total exon rows: {len(exon_bed_rows)}, intron rows: {len(intron_bed_rows)}")
