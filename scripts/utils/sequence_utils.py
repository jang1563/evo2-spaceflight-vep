"""
Sequence utilities for Evo2 VEP scoring.

Handles:
- Reference genome sequence extraction (via pyfaidx)
- Window construction around variants (centered)
- Variant sequence generation (SNVs and indels)
- Reference window deduplication for efficiency
- Reverse complement handling
"""

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from pyfaidx import Fasta
except ImportError:
    Fasta = None


# =============================================================================
# Constants
# =============================================================================

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
BASES = ["A", "C", "G", "T"]

# Default reference genome path — import from config if available
try:
    from .config import GENOME_PATH as _GENOME_PATH
    DEFAULT_GENOME = str(_GENOME_PATH)
except (ImportError, SystemError):
    try:
        from utils.config import GENOME_PATH as _GENOME_PATH
        DEFAULT_GENOME = str(_GENOME_PATH)
    except ImportError:
        DEFAULT_GENOME = "/path/to/your/scratch/evo2/data/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class Variant:
    """A single genomic variant."""
    chrom: str
    pos: int          # 0-based position
    ref: str
    alt: str
    variant_id: str = ""
    gene: str = ""
    region_type: str = ""  # "coding", "intronic", "promoter", "utr5", "utr3", "intergenic"
    clinvar_class: str = ""  # "Pathogenic", "Benign", "VUS", etc.
    clinvar_stars: int = 0
    dms_score: Optional[float] = None

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1

    @property
    def is_insertion(self) -> bool:
        return len(self.ref) < len(self.alt)

    @property
    def is_deletion(self) -> bool:
        return len(self.ref) > len(self.alt)

    @property
    def is_indel(self) -> bool:
        return self.is_insertion or self.is_deletion

    @property
    def is_frameshift(self) -> bool:
        """Check if indel causes frameshift (length not divisible by 3)."""
        if not self.is_indel:
            return False
        diff = abs(len(self.alt) - len(self.ref))
        return diff % 3 != 0

    @property
    def indel_length(self) -> int:
        """Signed indel length (positive=insertion, negative=deletion)."""
        return len(self.alt) - len(self.ref)

    def key(self) -> str:
        """Unique key for this variant."""
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass
class ScoringWindow:
    """A sequence window for Evo2 scoring."""
    chrom: str
    start: int
    end: int
    sequence: str
    variant: Optional[Variant] = None
    is_reference: bool = True
    window_hash: str = ""

    def __post_init__(self):
        if not self.window_hash:
            self.window_hash = hashlib.md5(self.sequence.encode()).hexdigest()[:12]


@dataclass
class ScoringResult:
    """Result of scoring a single variant."""
    variant: Variant
    ref_score: float
    alt_score: float
    delta: float = 0.0
    window_size: int = 0

    def __post_init__(self):
        self.delta = self.alt_score - self.ref_score


# =============================================================================
# Genome Access
# =============================================================================

class GenomeAccessor:
    """Thread-safe genome sequence accessor using pyfaidx."""

    def __init__(self, genome_path: str = DEFAULT_GENOME):
        if Fasta is None:
            raise ImportError("pyfaidx required: pip install pyfaidx")
        self.genome_path = genome_path
        self._fasta = None

    @property
    def fasta(self):
        if self._fasta is None:
            self._fasta = Fasta(self.genome_path)
        return self._fasta

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch sequence from genome. Coords are 0-based, half-open."""
        # pyfaidx uses 0-based indexing with slicing
        return str(self.fasta[chrom][start:end]).upper()

    def get_base(self, chrom: str, pos: int) -> str:
        """Get single base at position (0-based)."""
        return str(self.fasta[chrom][pos:pos + 1]).upper()

    def close(self):
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None


# =============================================================================
# Window Construction
# =============================================================================

def build_scoring_window(
    genome: GenomeAccessor,
    variant: Variant,
    window_size: int = 8192,
) -> Tuple[ScoringWindow, ScoringWindow]:
    """
    Build reference and alternate scoring windows centered on variant.

    For SNVs: simply substitute the base.
    For insertions: insert bases, trim downstream to maintain window_size.
    For deletions: remove bases, extend downstream to maintain window_size.

    Returns:
        (ref_window, alt_window)
    """
    half = window_size // 2
    center = variant.pos

    # Calculate window boundaries
    win_start = center - half
    win_end = center + half

    # Ensure we don't go below 0
    if win_start < 0:
        win_end += abs(win_start)
        win_start = 0

    # Fetch reference sequence (with extra buffer for indels)
    buffer = max(abs(variant.indel_length) + 100, 200) if variant.is_indel else 0
    ref_seq_extended = genome.fetch(variant.chrom, win_start, win_end + buffer)

    # Check we got enough sequence (pyfaidx silently truncates near chromosome ends)
    min_needed = window_size + (abs(variant.indel_length) if variant.is_indel else 0)
    if len(ref_seq_extended) < min_needed:
        raise ValueError(
            f"Cannot build {window_size}bp window for {variant.key()}: "
            f"only {len(ref_seq_extended)}bp available at "
            f"{variant.chrom}:{win_start}-{win_end + buffer} "
            f"(near chromosome boundary)"
        )

    # Position of variant within the window
    var_offset = center - win_start

    # Verify reference allele matches
    ref_at_pos = ref_seq_extended[var_offset:var_offset + len(variant.ref)]
    if ref_at_pos != variant.ref:
        raise ValueError(
            f"Reference mismatch at {variant.key()}: "
            f"expected '{variant.ref}', got '{ref_at_pos}'"
        )

    # Build reference window (exact window_size)
    ref_seq = ref_seq_extended[:window_size]

    # Build alternate window
    if variant.is_snv:
        alt_seq = ref_seq[:var_offset] + variant.alt + ref_seq[var_offset + 1:]
    elif variant.is_insertion:
        # Insert bases, then trim to window_size
        alt_extended = (
            ref_seq_extended[:var_offset + len(variant.ref)]
            + variant.alt[len(variant.ref):]
            + ref_seq_extended[var_offset + len(variant.ref):]
        )
        alt_seq = alt_extended[:window_size]
    elif variant.is_deletion:
        # Remove bases, then extend downstream to maintain window_size
        del_len = len(variant.ref) - len(variant.alt)
        alt_extended = (
            ref_seq_extended[:var_offset + len(variant.alt)]
            + ref_seq_extended[var_offset + len(variant.ref):]
        )
        alt_seq = alt_extended[:window_size]
    else:
        # MNV (multi-nucleotide variant): treat as block substitution
        alt_seq = (
            ref_seq[:var_offset]
            + variant.alt
            + ref_seq[var_offset + len(variant.ref):]
        )[:window_size]

    assert len(ref_seq) == window_size, f"ref_seq length {len(ref_seq)} != {window_size}"
    assert len(alt_seq) == window_size, f"alt_seq length {len(alt_seq)} != {window_size}"

    ref_window = ScoringWindow(
        chrom=variant.chrom, start=win_start, end=win_end,
        sequence=ref_seq, variant=variant, is_reference=True,
    )
    alt_window = ScoringWindow(
        chrom=variant.chrom, start=win_start, end=win_end,
        sequence=alt_seq, variant=variant, is_reference=False,
    )

    return ref_window, alt_window


# =============================================================================
# Variant Generation
# =============================================================================

def generate_all_snvs(
    genome: GenomeAccessor,
    chrom: str,
    start: int,
    end: int,
    gene: str = "",
) -> List[Variant]:
    """Generate all possible SNVs in a genomic region."""
    ref_seq = genome.fetch(chrom, start, end)
    variants = []
    for i, ref_base in enumerate(ref_seq):
        if ref_base not in BASES:
            continue
        pos = start + i
        for alt_base in BASES:
            if alt_base == ref_base:
                continue
            variants.append(Variant(
                chrom=chrom, pos=pos, ref=ref_base, alt=alt_base,
                variant_id=f"{chrom}:{pos}:{ref_base}>{alt_base}",
                gene=gene,
            ))
    return variants


def generate_radiation_snvs(
    genome: GenomeAccessor,
    chrom: str,
    start: int,
    end: int,
    gene: str = "",
) -> List[Variant]:
    """
    Generate radiation-characteristic SNVs (oxidative damage signature).

    Solar proton / trapped proton → C>A and G>T transversions via 8-oxoG (SBS18).
    NOT C>T (that's UV/SBS7).
    """
    ref_seq = genome.fetch(chrom, start, end)
    variants = []
    for i, ref_base in enumerate(ref_seq):
        pos = start + i
        if ref_base == "C":
            variants.append(Variant(
                chrom=chrom, pos=pos, ref="C", alt="A",
                variant_id=f"{chrom}:{pos}:C>A:radiation",
                gene=gene, region_type="radiation_C>A",
            ))
        elif ref_base == "G":
            variants.append(Variant(
                chrom=chrom, pos=pos, ref="G", alt="T",
                variant_id=f"{chrom}:{pos}:G>T:radiation",
                gene=gene, region_type="radiation_G>T",
            ))
    return variants


def generate_microhomology_deletions(
    genome: GenomeAccessor,
    chrom: str,
    start: int,
    end: int,
    gene: str = "",
    max_del_length: int = 5,
    max_mh_length: int = 4,
) -> List[Variant]:
    """
    Generate small deletions with microhomology (radiation indel signature).

    IR-induced deletions characteristically have 1-4bp microhomology at junctions.
    """
    ref_seq = genome.fetch(chrom, start, end + max_del_length + max_mh_length)
    variants = []

    for i in range(len(ref_seq) - max_del_length - max_mh_length):
        pos = start + i
        for del_len in range(1, max_del_length + 1):
            # Check for microhomology: bases after deletion match bases before
            deleted = ref_seq[i + 1:i + 1 + del_len]
            following = ref_seq[i + 1 + del_len:i + 1 + del_len + len(deleted)]

            if not following:
                continue

            # Count microhomology length
            mh_len = 0
            for a, b in zip(deleted, following):
                if a == b:
                    mh_len += 1
                else:
                    break

            if 1 <= mh_len <= max_mh_length:
                ref_allele = ref_seq[i:i + 1 + del_len]
                alt_allele = ref_seq[i:i + 1]
                if pos + len(ref_allele) > end:
                    continue
                variants.append(Variant(
                    chrom=chrom, pos=pos,
                    ref=ref_allele, alt=alt_allele,
                    variant_id=f"{chrom}:{pos}:del{del_len}:mh{mh_len}",
                    gene=gene, region_type=f"radiation_del_mh{mh_len}",
                ))

    return variants


# =============================================================================
# Reference Deduplication
# =============================================================================

def deduplicate_references(
    variants: List[Variant],
    window_size: int = 8192,
) -> Dict[str, List[Variant]]:
    """
    Group variants that share the same reference window.

    Variants at the same position with different alt alleles share
    one reference window. This saves ~35% compute for saturation mutagenesis.

    Returns:
        Dict mapping ref_window_key -> list of variants sharing that window.
    """
    groups: Dict[str, List[Variant]] = {}
    half = window_size // 2

    for var in variants:
        # Reference window is determined by position and window size
        win_start = max(0, var.pos - half)
        key = f"{var.chrom}:{win_start}:{win_start + window_size}"
        if key not in groups:
            groups[key] = []
        groups[key].append(var)

    return groups


# =============================================================================
# Annotation
# =============================================================================

def annotate_variant_region(variant: Variant, gene) -> str:
    """Annotate whether a variant falls in coding, intronic, promoter, or UTR region."""
    pos = variant.pos

    # Check if in any exon
    for exon in gene.exons:
        if exon.start <= pos < exon.end:
            # In exon — is it coding or UTR?
            if gene.cds_start <= pos < gene.cds_end:
                return "coding"
            elif (gene.strand == "+" and pos < gene.cds_start) or \
                 (gene.strand == "-" and pos >= gene.cds_end):
                return "utr5"
            else:
                return "utr3"

    # Check promoter
    if gene.promoter_start is not None and gene.promoter_end is not None:
        if gene.promoter_start <= pos < gene.promoter_end:
            return "promoter"

    # Otherwise intronic (within gene body) or intergenic
    if gene.start <= pos < gene.end:
        return "intronic"

    return "intergenic"


# =============================================================================
# Checkpointing
# =============================================================================

class ScoringCheckpoint:
    """Simple checkpoint manager for resumable scoring."""

    def __init__(self, checkpoint_dir: str):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.results_file = self.checkpoint_dir / "results.jsonl"
        self.scored_keys: set = set()
        self._load_existing()

    def _load_existing(self):
        """Load previously scored variant keys."""
        if self.results_file.exists():
            with open(self.results_file) as f:
                for line in f:
                    try:
                        data = json.loads(line)
                        self.scored_keys.add(data["variant_key"])
                    except (json.JSONDecodeError, KeyError):
                        continue

    def is_scored(self, variant: Variant) -> bool:
        return variant.key() in self.scored_keys

    def save_result(self, result: ScoringResult):
        """Append a scoring result."""
        data = {
            "variant_key": result.variant.key(),
            "gene": result.variant.gene,
            "chrom": result.variant.chrom,
            "pos": result.variant.pos,
            "ref": result.variant.ref,
            "alt": result.variant.alt,
            "region_type": result.variant.region_type,
            "ref_score": float(result.ref_score),
            "alt_score": float(result.alt_score),
            "delta": float(result.delta),
            "window_size": result.window_size,
            "clinvar_class": result.variant.clinvar_class,
            "clinvar_stars": result.variant.clinvar_stars,
            "dms_score": result.variant.dms_score,
        }
        with open(self.results_file, "a") as f:
            f.write(json.dumps(data) + "\n")
        self.scored_keys.add(result.variant.key())

    def load_results(self) -> List[dict]:
        """Load all results as list of dicts."""
        results = []
        if self.results_file.exists():
            with open(self.results_file) as f:
                for line in f:
                    try:
                        results.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
        return results

    @property
    def n_scored(self) -> int:
        return len(self.scored_keys)
