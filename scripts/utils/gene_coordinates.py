"""
Gene coordinates and exon boundaries for the 10-gene spaceflight VEP panel.

All coordinates are GRCh38/hg38.
Sources: NCBI RefSeq, Ensembl, UCSC Genome Browser (verified March 2026).

Gene panel:
  - 4 positive controls (with DMS data): BRCA1, TP53, CHEK2, DNMT3A
  - 6 novel targets: TERT, ATM, NFE2L2, CLOCK, MSTN, RAD51
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class Exon:
    """A single exon with genomic coordinates."""
    number: int
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Gene:
    """Gene with full genomic annotation."""
    symbol: str
    name: str
    chrom: str
    start: int          # 0-based, inclusive (gene body)
    end: int            # 0-based, exclusive (gene body)
    strand: str         # '+' or '-'
    refseq: str         # RefSeq transcript ID
    exons: List[Exon]   # sorted by genomic position
    cds_start: int = 0  # coding sequence start
    cds_end: int = 0    # coding sequence end

    # Metadata
    is_control: bool = False
    dms_source: str = ""
    spaceflight_evidence: str = ""
    rutter_classification: str = ""  # from Rutter/Mason 2024
    clinvar_vus_count: int = 0

    # Promoter region (for genes where it's important)
    promoter_start: Optional[int] = None
    promoter_end: Optional[int] = None

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def cds_length(self) -> int:
        """Total CDS length in bp (sum of coding exon portions)."""
        total = 0
        for exon in self.exons:
            ex_start = max(exon.start, self.cds_start)
            ex_end = min(exon.end, self.cds_end)
            if ex_start < ex_end:
                total += ex_end - ex_start
        return total

    @property
    def coding_exons(self) -> List[Exon]:
        """Exons that overlap the CDS."""
        return [e for e in self.exons if e.end > self.cds_start and e.start < self.cds_end]

    def get_region(self, padding: int = 0) -> Tuple[str, int, int]:
        """Get genomic region with optional padding."""
        return (self.chrom, max(0, self.start - padding), self.end + padding)

    def get_promoter_region(self, upstream: int = 2000, downstream: int = 500) -> Tuple[str, int, int]:
        """Get promoter region relative to TSS."""
        if self.promoter_start is not None and self.promoter_end is not None:
            return (self.chrom, self.promoter_start, self.promoter_end)
        if self.strand == '+':
            tss = self.start
            return (self.chrom, max(0, tss - upstream), tss + downstream)
        else:
            tss = self.end
            return (self.chrom, max(0, tss - downstream), tss + upstream)

    def region_str(self, padding: int = 0) -> str:
        """Format as chr:start-end string."""
        chrom, start, end = self.get_region(padding)
        return f"{chrom}:{start}-{end}"


# =============================================================================
# Gene Definitions
# =============================================================================
# Exon coordinates from NCBI RefSeq GFF3 (GCF_000001405.40_GRCh38.p14)
# Verified 2026-03-22 by parsing GFF3 annotation file.
# All coordinates are 0-based, half-open [start, end).
# Exons sorted by ascending genomic position.
# =============================================================================

TERT = Gene(
    symbol="TERT",
    name="Telomerase reverse transcriptase",
    chrom="chr5",
    start=1253166,
    end=1295068,
    strand="-",
    refseq="NM_198253.3",
    exons=[
        Exon(1,  1253166, 1253831),
        Exon(2,  1254367, 1254505),
        Exon(3,  1255286, 1255411),
        Exon(4,  1258597, 1258659),
        Exon(5,  1260473, 1260600),
        Exon(6,  1264403, 1264592),
        Exon(7,  1266463, 1266535),
        Exon(8,  1268519, 1268633),
        Exon(9,  1271118, 1271204),
        Exon(10, 1272184, 1272280),
        Exon(11, 1278640, 1278796),
        Exon(12, 1279290, 1279470),
        Exon(13, 1280157, 1280338),
        Exon(14, 1282428, 1282624),
        Exon(15, 1293312, 1294666),
        Exon(16, 1294770, 1295068),
    ],
    cds_start=1253727,
    cds_end=1294989,
    is_control=False,
    spaceflight_evidence="14.5% telomere elongation in NASA Twins Study; "
                         "core spaceflight biomarker",
    rutter_classification="Protective",
    clinvar_vus_count=180,
    # TERT promoter hotspots: C228T (chr5:1295113), C250T (chr5:1295135)
    promoter_start=1294568,  # TSS - 500bp
    promoter_end=1295568,    # TSS + 500bp (includes hotspots)
)

ATM = Gene(
    symbol="ATM",
    name="ATM serine/threonine kinase",
    chrom="chr11",
    start=108223066,
    end=108369102,
    strand="+",
    refseq="NM_000051.4",
    exons=[
        Exon(1,  108223066, 108223186),
        Exon(2,  108227594, 108227696),
        Exon(3,  108227775, 108227888),
        Exon(4,  108229177, 108229323),
        Exon(5,  108235669, 108235834),
        Exon(6,  108243952, 108244118),
        Exon(7,  108244787, 108245026),
        Exon(8,  108246963, 108247127),
        Exon(9,  108248932, 108249102),
        Exon(10, 108250700, 108251072),
        Exon(11, 108251836, 108252031),
        Exon(12, 108252816, 108252912),
        Exon(13, 108253813, 108254039),
        Exon(14, 108256214, 108256340),
        Exon(15, 108257480, 108257606),
        Exon(16, 108258985, 108259075),
        Exon(17, 108267170, 108267342),
        Exon(18, 108268409, 108268609),
        Exon(19, 108271063, 108271146),
        Exon(20, 108271250, 108271406),
        Exon(21, 108272531, 108272607),
        Exon(22, 108272721, 108272852),
        Exon(23, 108279490, 108279608),
        Exon(24, 108280994, 108281168),
        Exon(25, 108282709, 108282879),
        Exon(26, 108284226, 108284473),
        Exon(27, 108287599, 108287715),
        Exon(28, 108288976, 108289103),
        Exon(29, 108289601, 108289801),
        Exon(30, 108292618, 108292793),
        Exon(31, 108293312, 108293477),
        Exon(32, 108294926, 108295059),
        Exon(33, 108297286, 108297382),
        Exon(34, 108299713, 108299885),
        Exon(35, 108301647, 108301789),
        Exon(36, 108302852, 108303029),
        Exon(37, 108304674, 108304852),
        Exon(38, 108307896, 108307984),
        Exon(39, 108310159, 108310315),
        Exon(40, 108312410, 108312498),
        Exon(41, 108315822, 108315911),
        Exon(42, 108316010, 108316113),
        Exon(43, 108317372, 108317521),
        Exon(44, 108319953, 108320058),
        Exon(45, 108321300, 108321420),
        Exon(46, 108325309, 108325544),
        Exon(47, 108326057, 108326225),
        Exon(48, 108327644, 108327758),
        Exon(49, 108329020, 108329238),
        Exon(50, 108330213, 108330421),
        Exon(51, 108331443, 108331557),
        Exon(52, 108331878, 108332037),
        Exon(53, 108332761, 108332900),
        Exon(54, 108333885, 108333968),
        Exon(55, 108334968, 108335109),
        Exon(56, 108335844, 108335961),
        Exon(57, 108343221, 108343371),
        Exon(58, 108345742, 108345908),
        Exon(59, 108347278, 108347365),
        Exon(60, 108353765, 108353880),
        Exon(61, 108354810, 108354874),
        Exon(62, 108365081, 108365218),
        Exon(63, 108365324, 108369102),
    ],
    cds_start=108227624,
    cds_end=108365508,
    is_control=False,
    spaceflight_evidence="DSB sensing for GCR/HZE damage; key radiation response kinase",
    rutter_classification="Risk (LOF)",
    clinvar_vus_count=1500,
)

NFE2L2 = Gene(
    symbol="NFE2L2",
    name="NFE2 like bZIP transcription factor 2 (NRF2)",
    chrom="chr2",
    start=177230307,
    end=177264727,
    strand="-",
    refseq="NM_006164.5",
    exons=[
        Exon(1, 177230307, 177232008),
        Exon(2, 177232391, 177232583),
        Exon(3, 177233249, 177233339),
        Exon(4, 177234004, 177234271),
        Exon(5, 177264531, 177264727),
    ],
    cds_start=177230784,
    cds_end=177264576,
    is_control=False,
    spaceflight_evidence="4 ISS mouse papers: muscle, metabolism, immune, TMA; "
                         "master antioxidant TF",
    rutter_classification="Not listed",
    clinvar_vus_count=50,
)

CLOCK = Gene(
    symbol="CLOCK",
    name="Clock circadian regulator",
    chrom="chr4",
    start=55427902,
    end=55546909,
    strand="-",
    refseq="NM_004898.4",
    exons=[
        Exon(1,  55427902, 55435594),
        Exon(2,  55438281, 55438537),
        Exon(3,  55442431, 55442634),
        Exon(4,  55443686, 55443896),
        Exon(5,  55444632, 55444785),
        Exon(6,  55448778, 55448868),
        Exon(7,  55449395, 55449496),
        Exon(8,  55450090, 55450232),
        Exon(9,  55453053, 55453129),
        Exon(10, 55453676, 55453824),
        Exon(11, 55455896, 55456003),
        Exon(12, 55456217, 55456300),
        Exon(13, 55458891, 55459010),
        Exon(14, 55459147, 55459261),
        Exon(15, 55463684, 55463805),
        Exon(16, 55470716, 55470806),
        Exon(17, 55475962, 55476054),
        Exon(18, 55478814, 55478963),
        Exon(19, 55479639, 55479699),
        Exon(20, 55482738, 55482828),
        Exon(21, 55489373, 55489465),
        Exon(22, 55509911, 55510065),
        Exon(23, 55546781, 55546909),
    ],
    cds_start=55435414,
    cds_end=55482785,
    is_control=False,
    spaceflight_evidence="16 sunrises/day on ISS disrupts circadian rhythm; "
                         "I4/Twins Study pathway enrichment",
    rutter_classification="Not listed",
    clinvar_vus_count=60,
)

MSTN = Gene(
    symbol="MSTN",
    name="Myostatin",
    chrom="chr2",
    start=190055699,
    end=190062729,
    strand="+",
    refseq="NM_005259.3",
    exons=[
        Exon(1, 190055699, 190057638),
        Exon(2, 190060061, 190060435),
        Exon(3, 190062223, 190062729),
    ],
    cds_start=190057257,
    cds_end=190062596,
    is_control=False,
    spaceflight_evidence="LOF → muscle hypertrophy; protective against microgravity "
                         "muscle wasting",
    rutter_classification="Protective",
    clinvar_vus_count=30,
)

RAD51 = Gene(
    symbol="RAD51",
    name="RAD51 recombinase",
    chrom="chr15",
    start=40695173,
    end=40732340,
    strand="+",
    refseq="NM_002875.5",
    exons=[
        Exon(1,  40695173, 40695425),
        Exon(2,  40698756, 40698845),
        Exon(3,  40701063, 40701201),
        Exon(4,  40706176, 40706294),
        Exon(5,  40709024, 40709116),
        Exon(6,  40718804, 40718899),
        Exon(7,  40728710, 40728824),
        Exon(8,  40729504, 40729634),
        Exon(9,  40729852, 40729974),
        Exon(10, 40731054, 40732340),
    ],
    cds_start=40698758,
    cds_end=40731178,
    is_control=False,
    spaceflight_evidence="Core HR recombinase; strand invasion at DSBs from HZE; "
                         "loaded by BRCA2 (connects to BRCA1 control)",
    rutter_classification="Related (RAD50 in Rutter)",
    clinvar_vus_count=40,
)

BRCA1 = Gene(
    symbol="BRCA1",
    name="BRCA1 DNA repair associated",
    chrom="chr17",
    start=43044294,
    end=43125364,
    strand="-",
    refseq="NM_007294.4",
    exons=[
        Exon(1,  43044294, 43045802),
        Exon(2,  43047642, 43047703),
        Exon(3,  43049120, 43049194),
        Exon(4,  43051062, 43051117),
        Exon(5,  43057051, 43057135),
        Exon(6,  43063332, 43063373),
        Exon(7,  43063873, 43063951),
        Exon(8,  43067607, 43067695),
        Exon(9,  43070927, 43071238),
        Exon(10, 43074330, 43074521),
        Exon(11, 43076487, 43076614),
        Exon(12, 43082403, 43082575),
        Exon(13, 43090943, 43091032),
        Exon(14, 43091434, 43094860),  # Large exon (3,426 bp)
        Exon(15, 43095845, 43095922),
        Exon(16, 43097243, 43097289),
        Exon(17, 43099774, 43099880),
        Exon(18, 43104121, 43104261),
        Exon(19, 43104867, 43104956),
        Exon(20, 43106455, 43106533),
        Exon(21, 43115725, 43115779),
        Exon(22, 43124016, 43124115),
        Exon(23, 43125270, 43125364),
    ],
    cds_start=43045677,
    cds_end=43124096,
    is_control=True,
    dms_source="Findlay 2018 (SGE); MaveDB urn:mavedb:00000097-0-2; "
               "ClinGen PS3/BS3 approved. Coverage: exons 18-23 (RING+BRCT)",
    spaceflight_evidence="Rutter 2024 risk allele; DNA repair key for radiation response",
    rutter_classification="Risk",
    clinvar_vus_count=2100,
)

TP53 = Gene(
    symbol="TP53",
    name="Tumor protein p53",
    chrom="chr17",
    start=7668420,
    end=7687490,
    strand="-",
    refseq="NM_000546.6",
    exons=[
        Exon(1,  7668420, 7669690),
        Exon(2,  7670608, 7670715),
        Exon(3,  7673534, 7673608),
        Exon(4,  7673700, 7673837),
        Exon(5,  7674180, 7674290),
        Exon(6,  7674858, 7674971),
        Exon(7,  7675052, 7675236),
        Exon(8,  7675993, 7676272),
        Exon(9,  7676381, 7676403),
        Exon(10, 7676520, 7676622),
        Exon(11, 7687376, 7687490),
    ],
    cds_start=7669608,
    cds_end=7676594,
    is_control=True,
    dms_source="Giacomelli 2018 (nutlin-3 growth); MaveDB urn:mavedb:00000068-b-1; "
               "Full CDS coverage",
    spaceflight_evidence="Most mutated gene in astronaut CH (7 variants in 14 astronauts); "
                         "I4 crew mutations",
    rutter_classification="Risk",
    clinvar_vus_count=1200,
)

CHEK2 = Gene(
    symbol="CHEK2",
    name="Checkpoint kinase 2",
    chrom="chr22",
    start=28687742,
    end=28741820,
    strand="-",
    refseq="NM_007194.4",
    exons=[
        Exon(1,  28687742, 28687986),
        Exon(2,  28689134, 28689215),
        Exon(3,  28694031, 28694117),
        Exon(4,  28695126, 28695242),
        Exon(5,  28695709, 28695873),
        Exon(6,  28696900, 28696987),
        Exon(7,  28699837, 28699937),
        Exon(8,  28703504, 28703566),
        Exon(9,  28710005, 28710059),
        Exon(10, 28711908, 28712017),
        Exon(11, 28719394, 28719485),
        Exon(12, 28724976, 28725124),
        Exon(13, 28725242, 28725367),
        Exon(14, 28734402, 28734727),
        Exon(15, 28741768, 28741820),
    ],
    cds_start=28687896,
    cds_end=28734721,
    is_control=True,
    dms_source="McCarthy-Leo 2024 (4,885 SNVs); Kinase/SCD assay; Full CDS",
    spaceflight_evidence="DNA damage checkpoint kinase; activated by ATM",
    rutter_classification="Risk (inferred)",
    clinvar_vus_count=800,
)

DNMT3A = Gene(
    symbol="DNMT3A",
    name="DNA methyltransferase 3 alpha",
    chrom="chr2",
    start=25232960,
    end=25342590,
    strand="-",
    refseq="NM_175629.2",
    exons=[
        Exon(1,  25232960, 25234420),
        Exon(2,  25235706, 25235825),
        Exon(3,  25236935, 25237005),
        Exon(4,  25239129, 25239215),
        Exon(5,  25240301, 25240450),
        Exon(6,  25240639, 25240730),
        Exon(7,  25241561, 25241707),
        Exon(8,  25243897, 25243982),
        Exon(9,  25244154, 25244338),
        Exon(10, 25244539, 25244652),
        Exon(11, 25245252, 25245332),
        Exon(12, 25246019, 25246064),
        Exon(13, 25246159, 25246309),
        Exon(14, 25246619, 25246776),
        Exon(15, 25247050, 25247158),
        Exon(16, 25247590, 25247749),
        Exon(17, 25248036, 25248252),
        Exon(18, 25274940, 25275087),
        Exon(19, 25275499, 25275543),
        Exon(20, 25282440, 25282711),
        Exon(21, 25300138, 25300243),
        Exon(22, 25313912, 25314161),
        Exon(23, 25342429, 25342590),
    ],
    cds_start=25234278,
    cds_end=25313984,
    is_control=True,
    dms_source="Garcia et al. 2025 (paired DMS); bioRxiv 10.1101/2025.09.24.678339; "
               "2,036 variants across 3 libraries (alpha-2 helix, TRD, core RD interface)",
    spaceflight_evidence="6 mutations in 14 astronauts (2nd most after TP53); "
                         "found in I4 crew; clonal hematopoiesis marker",
    rutter_classification="Risk (CH)",
    clinvar_vus_count=375,
)


# =============================================================================
# Gene Registry
# =============================================================================

GENES: Dict[str, Gene] = {
    "TERT": TERT,
    "ATM": ATM,
    "NFE2L2": NFE2L2,
    "CLOCK": CLOCK,
    "MSTN": MSTN,
    "RAD51": RAD51,
    "BRCA1": BRCA1,
    "TP53": TP53,
    "CHEK2": CHEK2,
    "DNMT3A": DNMT3A,
}

CONTROL_GENES = [g for g in GENES.values() if g.is_control]
NOVEL_GENES = [g for g in GENES.values() if not g.is_control]

# TERT promoter hotspot positions (GRCh38, 0-based)
# These are on the minus strand (TERT coding strand); forward strand shows G
# C228T: minus strand C→T = forward strand G→A
# C250T: minus strand C→T = forward strand G→A
# NOTE: Verify exact hg38 positions against Kircher 2019 MPRA data in Phase 0.
# Positions below are from hg19 literature — may need coordinate adjustment.
TERT_PROMOTER_HOTSPOTS = {
    "C228T": ("chr5", 1295228, "G", "A"),  # Forward strand notation; creates de novo ETS site
    "C250T": ("chr5", 1295250, "G", "A"),  # Forward strand notation; creates de novo ETS site
}

# Astronaut-specific variants (from Brojakowska 2022, Commun Biol)
# Format: (gene, chrom, pos, ref, alt, astronaut_count, note)
ASTRONAUT_VARIANTS = [
    # TP53 — 7 variants in 14 astronauts (most mutated)
    # Exact positions from Brojakowska 2022 supplementary
    # To be filled in Phase 0 from supplementary data download

    # DNMT3A — 6 variants in 14 astronauts (2nd most mutated)
    # R882H is the most common clonal hematopoiesis hotspot
    # To be filled in Phase 0 from supplementary data download
]


# =============================================================================
# Helper Functions
# =============================================================================

def get_gene(symbol: str) -> Gene:
    """Get gene by symbol (case-insensitive)."""
    key = symbol.upper()
    if key not in GENES:
        raise ValueError(f"Unknown gene: {symbol}. Available: {list(GENES.keys())}")
    return GENES[key]


def get_all_regions(padding: int = 5000) -> Dict[str, Tuple[str, int, int]]:
    """Get all gene regions with padding for sequence extraction."""
    return {sym: gene.get_region(padding) for sym, gene in GENES.items()}


def get_bed_entries(padding: int = 0) -> List[str]:
    """Generate BED-format entries for all genes."""
    entries = []
    for sym, gene in sorted(GENES.items(), key=lambda x: (x[1].chrom, x[1].start)):
        chrom, start, end = gene.get_region(padding)
        entries.append(f"{chrom}\t{start}\t{end}\t{sym}\t0\t{gene.strand}")
    return entries


def write_bed(filepath: str, padding: int = 5000) -> None:
    """Write BED file for all gene regions."""
    entries = get_bed_entries(padding)
    with open(filepath, 'w') as f:
        for entry in entries:
            f.write(entry + '\n')


def summary_table() -> str:
    """Print a summary table of all genes."""
    lines = [
        f"{'Gene':<10} {'Chr':<6} {'Start':>12} {'End':>12} {'Span':>8} "
        f"{'Strand':<3} {'Exons':>5} {'Control':<8} {'VUS':>5}"
    ]
    lines.append("-" * 80)
    for sym in ["BRCA1", "TP53", "CHEK2", "DNMT3A",
                "TERT", "ATM", "NFE2L2", "CLOCK", "MSTN", "RAD51"]:
        g = GENES[sym]
        span = f"{g.length / 1000:.1f}kb"
        ctrl = "Yes" if g.is_control else ""
        lines.append(
            f"{g.symbol:<10} {g.chrom:<6} {g.start:>12,} {g.end:>12,} {span:>8} "
            f"{g.strand:<3} {g.n_exons:>5} {ctrl:<8} {g.clinvar_vus_count:>5}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print("Evo2 Spaceflight VEP — Gene Panel")
    print("=" * 80)
    print(summary_table())
    print()
    print(f"Total genes: {len(GENES)}")
    print(f"  Controls (with DMS): {len(CONTROL_GENES)}")
    print(f"  Novel targets: {len(NOVEL_GENES)}")
    print()
    for g in CONTROL_GENES:
        print(f"  {g.symbol}: {g.dms_source[:60]}...")
