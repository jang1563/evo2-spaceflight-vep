"""
HGVS-to-genomic coordinate mapping for DMS data.

Maps DMS variant notations (HGVS coding/protein/non-coding) to
genomic coordinates (chrom, pos, ref, alt) for Evo2 scoring.

Supports:
  - c.XN>M: CDS-relative nucleotide substitutions
  - p.Xxx###Yyy: Protein-level amino acid substitutions (back-translation)
  - n.XN>M: Non-coding relative positions (e.g., TERT MPRA)
"""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .gene_coordinates import Gene

# Complement map for strand conversion
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Standard genetic code (codon -> amino acid)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse: amino acid -> list of codons
AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# Three-letter to one-letter amino acid codes
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}
AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()}


@dataclass
class CdsMap:
    """
    Maps CDS positions to genomic coordinates for a gene.

    CDS position 1 = first nucleotide of start codon.
    For minus-strand genes, genomic coordinates decrease as CDS positions increase.
    """
    gene: Gene
    # List of (cds_start_pos, cds_end_pos, genomic_start, genomic_end, exon_number)
    # cds_start_pos: 1-based CDS position of the first base in this exon segment
    # For plus-strand: genomic_start < genomic_end
    # For minus-strand: genomic_start > genomic_end (stored as the highest genomic pos first)
    segments: List[Tuple[int, int, int, int, int]]  # (cds_start, cds_end, geno_start, geno_end, exon_num)
    total_cds_length: int

    def cds_to_genomic(self, cds_pos: int) -> Tuple[str, int]:
        """
        Convert 1-based CDS position to 0-based genomic coordinate.

        Returns (chrom, genomic_pos) where genomic_pos is 0-based.
        Supports negative positions (5'UTR: c.-1, c.-2, ...).
        """
        # Handle 5'UTR positions (c.-1, c.-2, etc.)
        if cds_pos < 0:
            # c.-N is N bases upstream of c.1 (in the 5'UTR)
            N = abs(cds_pos)
            first_seg = self.segments[0]
            geno_start_c1 = first_seg[2]  # genomic pos of c.1
            if self.gene.strand == "+":
                return (self.gene.chrom, geno_start_c1 - N)
            else:
                return (self.gene.chrom, geno_start_c1 + N)

        if cds_pos < 1 or cds_pos > self.total_cds_length:
            raise ValueError(
                f"CDS position {cds_pos} out of range [1, {self.total_cds_length}] "
                f"for {self.gene.symbol}"
            )

        for seg_cds_start, seg_cds_end, geno_start, geno_end, _ in self.segments:
            if seg_cds_start <= cds_pos <= seg_cds_end:
                offset = cds_pos - seg_cds_start
                if self.gene.strand == "+":
                    return (self.gene.chrom, geno_start + offset)
                else:
                    # Minus strand: CDS goes in opposite direction
                    return (self.gene.chrom, geno_start - offset)

        raise ValueError(f"CDS position {cds_pos} not found in exon segments for {self.gene.symbol}")

    def cds_to_genomic_intronic(self, cds_pos: int, intronic_offset: int) -> Tuple[str, int]:
        """
        Convert CDS boundary position + intronic offset to 0-based genomic coordinate.

        For donor variants (c.XXX+N): cds_pos is the last CDS base of an exon.
        For acceptor variants (c.XXX-N): cds_pos is the first CDS base of an exon.

        Intronic bases use coding-strand orientation in HGVS notation.
        """
        N = abs(intronic_offset)

        if intronic_offset > 0:
            # Donor: cds_pos should be seg_cds_end (last base of an exon)
            for seg_cds_start, seg_cds_end, geno_start, geno_end, _ in self.segments:
                if cds_pos == seg_cds_end:
                    # geno_end = last genomic position of segment
                    # For +strand: geno_end is highest pos; intron goes higher
                    # For -strand: geno_end is lowest pos; intron goes lower
                    if self.gene.strand == "+":
                        return (self.gene.chrom, geno_end + N)
                    else:
                        return (self.gene.chrom, geno_end - N)
            raise ValueError(
                f"CDS position {cds_pos} is not at a donor splice boundary "
                f"for {self.gene.symbol}"
            )

        elif intronic_offset < 0:
            # Acceptor: cds_pos should be seg_cds_start (first base of an exon)
            for seg_cds_start, seg_cds_end, geno_start, geno_end, _ in self.segments:
                if cds_pos == seg_cds_start:
                    # geno_start = first genomic position of segment
                    # For +strand: geno_start is lowest pos; intron goes lower
                    # For -strand: geno_start is highest pos; intron goes higher
                    if self.gene.strand == "+":
                        return (self.gene.chrom, geno_start - N)
                    else:
                        return (self.gene.chrom, geno_start + N)
            raise ValueError(
                f"CDS position {cds_pos} is not at an acceptor splice boundary "
                f"for {self.gene.symbol}"
            )

        else:
            return self.cds_to_genomic(cds_pos)

    def codon_genomic_positions(self, aa_pos: int) -> List[Tuple[str, int]]:
        """
        Get the 3 genomic positions of a codon at amino acid position aa_pos (1-based).

        Returns list of 3 (chrom, genomic_pos) tuples.
        """
        cds_start = (aa_pos - 1) * 3 + 1  # 1-based CDS position of codon start
        return [self.cds_to_genomic(cds_start + i) for i in range(3)]


def build_cds_map(gene: Gene) -> CdsMap:
    """
    Build a CDS-to-genomic coordinate map for a gene.

    Handles both plus and minus strand genes correctly.
    """
    # Get coding exons and their CDS overlaps
    coding_parts = []
    for exon in gene.exons:
        # Overlap of this exon with the CDS region
        overlap_start = max(exon.start, gene.cds_start)
        overlap_end = min(exon.end, gene.cds_end)
        if overlap_start < overlap_end:
            coding_parts.append((overlap_start, overlap_end, exon.number))

    if not coding_parts:
        raise ValueError(f"No coding exons found for {gene.symbol}")

    # Sort by CDS order
    if gene.strand == "+":
        # Plus strand: CDS goes left to right (ascending genomic coords)
        coding_parts.sort(key=lambda x: x[0])
    else:
        # Minus strand: CDS goes right to left (descending genomic coords)
        coding_parts.sort(key=lambda x: x[0], reverse=True)

    # Build segments with CDS positions
    segments = []
    cds_pos = 1
    for geno_start, geno_end, exon_num in coding_parts:
        seg_length = geno_end - geno_start

        if gene.strand == "+":
            segments.append((cds_pos, cds_pos + seg_length - 1, geno_start, geno_end - 1, exon_num))
        else:
            # For minus strand, CDS starts at the highest genomic position
            segments.append((cds_pos, cds_pos + seg_length - 1, geno_end - 1, geno_start, exon_num))

        cds_pos += seg_length

    total_length = cds_pos - 1
    return CdsMap(gene=gene, segments=segments, total_cds_length=total_length)


def complement(base: str) -> str:
    """Get the complement of a DNA base."""
    return COMPLEMENT.get(base.upper(), base)


def reverse_complement(seq: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    return "".join(complement(b) for b in reversed(seq))


# =============================================================================
# HGVS Parsing
# =============================================================================

@dataclass
class ParsedHGVS:
    """Parsed HGVS variant notation."""
    type: str        # "c" (coding), "p" (protein), "n" (non-coding)
    position: int    # Position (1-based)
    ref: str         # Reference allele/amino acid
    alt: str         # Alternate allele/amino acid
    is_synonymous: bool = False
    raw: str = ""    # Original HGVS string
    intronic_offset: int = 0  # For c.XXX+N (positive) or c.XXX-N (negative)


def parse_hgvs_coding(hgvs: str) -> Optional[ParsedHGVS]:
    """
    Parse HGVS coding DNA notation like:
      c.38T>C              — exonic SNV
      c.80+3T>C            — intronic donor SNV
      c.81-5G>A            — intronic acceptor SNV
      NM_007294.3:c.38T>C  — with transcript prefix (stripped)

    Returns None for complex variants (indels, multi-nucleotide).
    Only handles simple SNVs (exonic and intronic).
    """
    hgvs = hgvs.strip()

    # Strip transcript prefix (e.g., "NM_007294.3:")
    m_prefix = re.match(r"[A-Z]{2}_\d+\.\d+:", hgvs)
    if m_prefix:
        hgvs = hgvs[m_prefix.end():]

    if not hgvs.startswith("c."):
        return None

    # Simple exonic SNV: c.38T>C
    m = re.match(r"c\.(\d+)([ACGT])>([ACGT])$", hgvs)
    if m:
        pos = int(m.group(1))
        ref = m.group(2)
        alt = m.group(3)
        return ParsedHGVS(type="c", position=pos, ref=ref, alt=alt, raw=hgvs)

    # Intronic SNV: c.80+3T>C (donor) or c.81-5G>A (acceptor)
    m = re.match(r"c\.(\d+)([+-])(\d+)([ACGT])>([ACGT])$", hgvs)
    if m:
        pos = int(m.group(1))
        sign = m.group(2)
        offset = int(m.group(3))
        ref = m.group(4)
        alt = m.group(5)
        intronic_offset = offset if sign == "+" else -offset
        return ParsedHGVS(
            type="c", position=pos, ref=ref, alt=alt,
            raw=hgvs, intronic_offset=intronic_offset,
        )

    # 5'UTR SNV: c.-1A>T (negative CDS position = upstream of start codon)
    m = re.match(r"c\.(-\d+)([ACGT])>([ACGT])$", hgvs)
    if m:
        pos = int(m.group(1))  # Negative integer
        ref = m.group(2)
        alt = m.group(3)
        return ParsedHGVS(type="c", position=pos, ref=ref, alt=alt, raw=hgvs)

    return None


def parse_hgvs_protein(hgvs: str) -> Optional[ParsedHGVS]:
    """
    Parse HGVS protein notation like:
      p.Val13Ala
      p.Asp391Glu
      p.Pro390Ter  (nonsense)
      p.Asp391=    (synonymous)

    Returns None for complex variants.
    """
    hgvs = hgvs.strip()
    if not hgvs.startswith("p."):
        return None

    # Synonymous: p.Asp391=
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)=$", hgvs)
    if m:
        aa3 = m.group(1)
        pos = int(m.group(2))
        aa1 = AA3_TO_1.get(aa3, "?")
        return ParsedHGVS(
            type="p", position=pos, ref=aa1, alt=aa1,
            is_synonymous=True, raw=hgvs
        )

    # Missense/nonsense: p.Val13Ala or p.Pro390Ter
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", hgvs)
    if m:
        ref_aa3 = m.group(1)
        pos = int(m.group(2))
        alt_aa3 = m.group(3)
        ref_aa1 = AA3_TO_1.get(ref_aa3, "?")
        alt_aa1 = AA3_TO_1.get(alt_aa3, "?")
        return ParsedHGVS(
            type="p", position=pos, ref=ref_aa1, alt=alt_aa1,
            is_synonymous=False, raw=hgvs
        )

    return None


def parse_hgvs_noncoding(hgvs: str) -> Optional[ParsedHGVS]:
    """
    Parse HGVS non-coding notation like:
      n.3G>A
      n.292C>A
      n.1=   (reference/wildtype)

    Returns None for complex variants.
    """
    hgvs = hgvs.strip()
    if not hgvs.startswith("n."):
        return None

    # Reference: n.1=
    m = re.match(r"n\.(\d+)=$", hgvs)
    if m:
        pos = int(m.group(1))
        return ParsedHGVS(type="n", position=pos, ref="=", alt="=",
                         is_synonymous=True, raw=hgvs)

    # SNV: n.3G>A
    m = re.match(r"n\.(\d+)([ACGT])>([ACGT])$", hgvs)
    if m:
        pos = int(m.group(1))
        ref = m.group(2)
        alt = m.group(3)
        return ParsedHGVS(type="n", position=pos, ref=ref, alt=alt, raw=hgvs)

    return None


# =============================================================================
# HGVS to Genomic Coordinate Conversion
# =============================================================================

@dataclass
class GenomicVariant:
    """A variant in genomic coordinates, ready for Evo2 scoring."""
    chrom: str
    pos: int         # 0-based genomic position
    ref: str         # Reference allele (forward strand)
    alt: str         # Alternate allele (forward strand)
    gene: str
    hgvs: str        # Original HGVS notation
    dms_score: float  # Original DMS functional score


def coding_hgvs_to_genomic(
    parsed: ParsedHGVS,
    cds_map: CdsMap,
) -> Optional[GenomicVariant]:
    """
    Convert a parsed coding HGVS variant to genomic coordinates.

    Handles both exonic (c.38T>C) and intronic (c.80+3T>C) variants.
    The HGVS ref/alt are in mRNA (coding strand) orientation.
    For minus-strand genes, we complement them to get forward-strand alleles.
    """
    try:
        if parsed.intronic_offset != 0:
            chrom, geno_pos = cds_map.cds_to_genomic_intronic(
                parsed.position, parsed.intronic_offset)
        else:
            chrom, geno_pos = cds_map.cds_to_genomic(parsed.position)
    except ValueError:
        return None

    gene = cds_map.gene

    # Convert mRNA alleles to forward-strand alleles
    if gene.strand == "-":
        ref_fwd = complement(parsed.ref)
        alt_fwd = complement(parsed.alt)
    else:
        ref_fwd = parsed.ref
        alt_fwd = parsed.alt

    return GenomicVariant(
        chrom=chrom,
        pos=geno_pos,
        ref=ref_fwd,
        alt=alt_fwd,
        gene=gene.symbol,
        hgvs=parsed.raw,
        dms_score=0.0,  # Filled in by caller
    )


def protein_hgvs_to_genomic_candidates(
    parsed: ParsedHGVS,
    cds_map: CdsMap,
    genome_accessor=None,
) -> List[GenomicVariant]:
    """
    Convert a protein-level HGVS variant to candidate genomic SNVs.

    Since one amino acid change can result from multiple nucleotide changes,
    this returns ALL possible single-nucleotide variants that produce the
    specified amino acid substitution.

    If genome_accessor is provided, filters to variants whose ref matches
    the actual reference genome.
    """
    gene = cds_map.gene

    # Get the 3 genomic positions of the codon
    try:
        codon_positions = cds_map.codon_genomic_positions(parsed.position)
    except ValueError:
        return []

    # For each position in the codon, try single-nucleotide changes
    # that convert ref_aa to alt_aa
    candidates = []

    # Get all codons for ref and alt amino acids
    ref_codons = AA_TO_CODONS.get(parsed.ref, [])
    alt_codons = AA_TO_CODONS.get(parsed.alt, [])

    if not ref_codons or not alt_codons:
        return []

    # For each possible ref codon that matches the reference amino acid,
    # find single-nt changes that produce any alt codon for the alt amino acid
    for ref_codon in ref_codons:
        for alt_codon in alt_codons:
            if ref_codon == alt_codon:
                continue
            # Count differences
            diffs = [(i, ref_codon[i], alt_codon[i])
                     for i in range(3) if ref_codon[i] != alt_codon[i]]
            # Only single-nucleotide changes
            if len(diffs) != 1:
                continue

            codon_idx, ref_base, alt_base = diffs[0]
            chrom, geno_pos = codon_positions[codon_idx]

            # Convert to forward strand
            if gene.strand == "-":
                ref_fwd = complement(ref_base)
                alt_fwd = complement(alt_base)
            else:
                ref_fwd = ref_base
                alt_fwd = alt_base

            candidates.append(GenomicVariant(
                chrom=chrom,
                pos=geno_pos,
                ref=ref_fwd,
                alt=alt_fwd,
                gene=gene.symbol,
                hgvs=parsed.raw,
                dms_score=0.0,
            ))

    # If we have a genome accessor, filter to candidates whose ref matches
    if genome_accessor is not None:
        verified = []
        for cand in candidates:
            actual_ref = genome_accessor.fetch(cand.chrom, cand.pos, cand.pos + 1)
            if actual_ref.upper() == cand.ref.upper():
                verified.append(cand)
        return verified

    return candidates


def noncoding_hgvs_to_genomic(
    parsed: ParsedHGVS,
    reference_start: int,
    chrom: str,
    strand: str,
    gene_symbol: str,
) -> Optional[GenomicVariant]:
    """
    Convert a non-coding HGVS variant to genomic coordinates.

    Args:
        parsed: Parsed HGVS notation (n.XN>M)
        reference_start: 0-based genomic start of the non-coding reference region
        chrom: Chromosome
        strand: '+' or '-'
        gene_symbol: Gene name for annotation

    For the TERT MPRA (Kircher 2019), the reference is a specific promoter region.
    Position n.1 maps to the first base of that region.
    """
    if parsed.is_synonymous:
        return None

    # Map position to genomic coordinate
    if strand == "+":
        geno_pos = reference_start + parsed.position - 1
        ref_fwd = parsed.ref
        alt_fwd = parsed.alt
    else:
        # Minus strand: positions increase as genomic coords decrease
        geno_pos = reference_start - parsed.position + 1
        ref_fwd = complement(parsed.ref)
        alt_fwd = complement(parsed.alt)

    return GenomicVariant(
        chrom=chrom,
        pos=geno_pos,
        ref=ref_fwd,
        alt=alt_fwd,
        gene=gene_symbol,
        hgvs=parsed.raw,
        dms_score=0.0,
    )


# =============================================================================
# Batch Processing for DMS Datasets
# =============================================================================

def map_brca1_dms(dms_rows: list, gene: Gene) -> List[GenomicVariant]:
    """
    Map BRCA1 Findlay 2018 SGE DMS data to genomic variants.

    Input: rows from MaveDB urn:mavedb:00000097-0-2 CSV with columns:
      - hgvs_nt: NM_007294.3:c.XXXN>M (with transcript prefix)
      - score: SGE functional score
    Handles both exonic SNVs and intronic splice-region variants.
    Transcript prefix is stripped automatically by parse_hgvs_coding().
    """
    cds_map = build_cds_map(gene)
    variants = []
    skipped = {"mnv": 0, "no_nt": 0, "no_score": 0, "parse_fail": 0, "map_fail": 0}

    for row in dms_rows:
        hgvs_nt = row.get("hgvs_nt", "")
        score_str = row.get("score", "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not hgvs_nt or hgvs_nt == "NA":
            skipped["no_nt"] += 1
            continue
        if ";" in hgvs_nt:
            skipped["mnv"] += 1
            continue

        parsed = parse_hgvs_coding(hgvs_nt)
        if parsed is None:
            skipped["parse_fail"] += 1
            continue

        gv = coding_hgvs_to_genomic(parsed, cds_map)
        if gv is None:
            skipped["map_fail"] += 1
            continue

        gv.dms_score = float(score_str)
        variants.append(gv)

    return variants


def map_tp53_dms(dms_rows: list, gene: Gene, genome_accessor=None) -> List[GenomicVariant]:
    """
    Map TP53 Giacomelli 2018 DMS data to genomic variants.

    Input: rows from CSV with columns: hgvs_pro, score
    Only has protein-level HGVS (no nucleotide notation).
    Returns multiple candidates per amino acid change (all possible SNVs).
    If genome_accessor is provided, filters to reference-matching candidates.
    """
    cds_map = build_cds_map(gene)
    variants = []
    skipped = {"no_pro": 0, "no_score": 0, "parse_fail": 0, "no_candidates": 0, "synonymous": 0}

    for row in dms_rows:
        hgvs_pro = row.get("hgvs_pro", "")
        score_str = row.get("score", "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not hgvs_pro or hgvs_pro == "NA":
            skipped["no_pro"] += 1
            continue

        parsed = parse_hgvs_protein(hgvs_pro)
        if parsed is None:
            skipped["parse_fail"] += 1
            continue
        if parsed.is_synonymous:
            skipped["synonymous"] += 1
            # Still map synonymous for calibration (should have delta ~0)
            # But we need to handle them differently since ref==alt in protein
            continue

        candidates = protein_hgvs_to_genomic_candidates(
            parsed, cds_map, genome_accessor=genome_accessor
        )
        if not candidates:
            skipped["no_candidates"] += 1
            continue

        for gv in candidates:
            gv.dms_score = float(score_str)
            variants.append(gv)

    return variants


def map_chek2_dms(dms_rows: list, gene: Gene, score_column: str = "RCS_this_SNV") -> List[GenomicVariant]:
    """
    Map CHEK2 McCarthy-Leo 2024 DMS data to genomic variants.

    Input: rows from CSV with columns: hgvs_nt, hgvs_pro, score, RCS_this_SNV
    Has nucleotide HGVS for direct mapping.

    Args:
        score_column: Column to use as DMS score. "RCS_this_SNV" is the continuous
                     Rad53 complementation score. "score" is binary (0/2).
    """
    cds_map = build_cds_map(gene)
    variants = []
    skipped = {"no_nt": 0, "no_score": 0, "parse_fail": 0, "map_fail": 0}

    for row in dms_rows:
        hgvs_nt = row.get("hgvs_nt", "")
        score_str = row.get(score_column, "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not hgvs_nt or hgvs_nt == "NA":
            skipped["no_nt"] += 1
            continue

        parsed = parse_hgvs_coding(hgvs_nt)
        if parsed is None:
            skipped["parse_fail"] += 1
            continue

        gv = coding_hgvs_to_genomic(parsed, cds_map)
        if gv is None:
            skipped["map_fail"] += 1
            continue

        gv.dms_score = float(score_str)
        variants.append(gv)

    return variants


def map_dnmt3a_dms(dms_rows: list, gene: Gene) -> List[GenomicVariant]:
    """
    Map DNMT3A stability assay data (Huang 2022 fallback; Garcia et al. 2025 preferred) to genomic variants.

    Input: rows from CSV with columns: Mutation (e.g., "R19W"),
           Stability ratio (Normalized to DNMT3AWT)
    Mutation format: single-letter ref AA + position + single-letter alt AA.
    Score: stability ratio (lower = more unstable/damaging).
    """
    cds_map = build_cds_map(gene)
    variants = []
    skipped = {"parse_fail": 0, "no_candidates": 0, "no_score": 0}

    for row in dms_rows:
        mutation = row.get("Mutation", "").strip()
        score_str = row.get("Stability ratio (Normalized to DNMT3AWT)", "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not mutation:
            skipped["parse_fail"] += 1
            continue

        # Parse mutation format: R19W → ref=R, pos=19, alt=W
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", mutation)
        if not m:
            skipped["parse_fail"] += 1
            continue

        ref_aa = m.group(1)
        aa_pos = int(m.group(2))
        alt_aa = m.group(3)

        if ref_aa == alt_aa:
            continue

        # Build a protein-level parsed HGVS
        parsed = ParsedHGVS(
            type="p", position=aa_pos, ref=ref_aa, alt=alt_aa,
            is_synonymous=False, raw=f"p.{ref_aa}{aa_pos}{alt_aa}",
        )

        candidates = protein_hgvs_to_genomic_candidates(parsed, cds_map)
        if not candidates:
            skipped["no_candidates"] += 1
            continue

        for gv in candidates:
            gv.dms_score = float(score_str)
            variants.append(gv)

    return variants


def map_dnmt3a_garcia_dms(
    dms_rows: list, gene: Gene, score_column: str = "score",
    genome_accessor=None,
) -> List[GenomicVariant]:
    """
    Map DNMT3A Garcia et al. 2025 methylation activity DMS data.

    Input: rows from processed CSV (dnmt3a_garcia_2025_wt.csv) with columns:
      - hgvs_pro: p.R4L format (1-letter AA codes)
      - pos: amino acid position
      - ref_aa, alt_aa: single-letter codes
      - score: methylation activity enrichment (0 = WT, negative = LOF)
      - stop_normalized_score: stop codon-normalized score
      - library: lib1/lib2/lib3

    Reference: Garcia, Lavidor et al. 2025 (bioRxiv 10.1101/2025.09.24.678339)
    """
    cds_map = build_cds_map(gene)
    variants = []
    skipped = {"parse_fail": 0, "no_candidates": 0, "no_score": 0}

    for row in dms_rows:
        ref_aa = row.get("ref_aa", "").strip()
        alt_aa = row.get("alt_aa", "").strip()
        pos_str = row.get("pos", "")
        score_str = row.get(score_column, "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not ref_aa or not alt_aa or not pos_str:
            skipped["parse_fail"] += 1
            continue

        aa_pos = int(pos_str)
        if ref_aa == alt_aa:
            continue

        parsed = ParsedHGVS(
            type="p", position=aa_pos, ref=ref_aa, alt=alt_aa,
            is_synonymous=False, raw=f"p.{ref_aa}{aa_pos}{alt_aa}",
        )

        candidates = protein_hgvs_to_genomic_candidates(
            parsed, cds_map, genome_accessor=genome_accessor
        )
        if not candidates:
            skipped["no_candidates"] += 1
            continue

        for gv in candidates:
            gv.dms_score = float(score_str)
            variants.append(gv)

    return variants


def map_tert_mpra(mpra_rows: list, reference_start: int, gene: Gene) -> List[GenomicVariant]:
    """
    Map TERT MPRA (Kircher 2019) data to genomic variants.

    The MPRA construct uses the FORWARD strand sequence, even though TERT
    is on the minus strand. The n. positions run 5'→3' on the forward strand:
    n.1 maps to reference_start, n.2 to reference_start+1, etc.
    Bases in the MPRA data are already in forward-strand orientation.

    Args:
        mpra_rows: Rows from CSV with columns: hgvs_nt, score
        reference_start: 0-based genomic position of n.1 (forward strand)
        gene: TERT gene object (used only for chrom and symbol)
    """
    variants = []
    skipped = {"no_nt": 0, "no_score": 0, "parse_fail": 0, "synonymous": 0}

    for row in mpra_rows:
        hgvs_nt = row.get("hgvs_nt", "")
        score_str = row.get("score", "")

        if not score_str or score_str == "NA":
            skipped["no_score"] += 1
            continue
        if not hgvs_nt or hgvs_nt == "NA":
            skipped["no_nt"] += 1
            continue

        parsed = parse_hgvs_noncoding(hgvs_nt)
        if parsed is None:
            skipped["parse_fail"] += 1
            continue
        if parsed.is_synonymous:
            skipped["synonymous"] += 1
            continue

        # Use strand="+" because the MPRA construct is on the forward strand
        gv = noncoding_hgvs_to_genomic(
            parsed, reference_start,
            chrom=gene.chrom, strand="+",
            gene_symbol=gene.symbol,
        )
        if gv is None:
            continue

        gv.dms_score = float(score_str)
        variants.append(gv)

    return variants
