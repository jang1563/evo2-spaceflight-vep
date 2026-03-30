#!/bin/bash
# =============================================================================
# Download all reference data for Evo2 Spaceflight VEP project
# Run on Cayuga HPC login node (has internet access)
#
# Usage: bash download_data.sh
# =============================================================================

set -euo pipefail

BASE="/athena/masonlab/scratch/users/jak4013/evo2/data"
PROJECT_DATA="/athena/masonlab/scratch/users/jak4013/evo2/project_spaceflight_vep/data"
TMPDIR="/athena/masonlab/scratch/users/jak4013/evo2/tmp"
export TMPDIR

echo "=== Evo2 VEP Data Download: $(date) ==="
echo "Base directory: $BASE"

# -------------------------------------------------------------------------
# 1. GRCh38 Reference Genome
# -------------------------------------------------------------------------
echo ""
echo "=== 1. GRCh38 Reference Genome ==="
mkdir -p "$BASE/GRCh38"
cd "$BASE/GRCh38"

if [ ! -f "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" ]; then
    echo "Downloading GRCh38 no-alt analysis set..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    echo "Decompressing..."
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    echo "Building faidx index..."
    # Need samtools or pyfaidx to index
    python3 -c "from pyfaidx import Fasta; f = Fasta('GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'); print(f'Indexed: {len(f.keys())} contigs')"
else
    echo "GRCh38 genome already exists, skipping."
fi

# -------------------------------------------------------------------------
# 2. ClinVar VCF (GRCh38)
# -------------------------------------------------------------------------
echo ""
echo "=== 2. ClinVar VCF ==="
mkdir -p "$BASE/clinvar"
cd "$BASE/clinvar"

echo "Downloading latest ClinVar VCF (GRCh38)..."
wget -q -O clinvar.vcf.gz "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
wget -q -O clinvar.vcf.gz.tbi "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
echo "ClinVar downloaded."

# Count variants for our genes
echo "Variant counts per gene (approximate):"
for gene in TERT ATM NFE2L2 CLOCK MSTN RAD51 BRCA1 TP53 CHEK2 DNMT3A; do
    count=$(zgrep -c "GENEINFO=$gene:" clinvar.vcf.gz 2>/dev/null || echo "0")
    echo "  $gene: $count"
done

# -------------------------------------------------------------------------
# 3. NCBI Gene Annotation (GFF3 for exact exon coordinates)
# -------------------------------------------------------------------------
echo ""
echo "=== 3. NCBI RefSeq GFF3 ==="
mkdir -p "$BASE/annotation"
cd "$BASE/annotation"

if [ ! -f "GCF_000001405.40_GRCh38.p14_genomic.gff.gz" ]; then
    echo "Downloading RefSeq GFF3 annotation..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
else
    echo "RefSeq GFF3 already exists, skipping."
fi

# -------------------------------------------------------------------------
# 4. DMS Datasets (MaveDB)
# -------------------------------------------------------------------------
echo ""
echo "=== 4. DMS Datasets (project-specific) ==="
mkdir -p "$PROJECT_DATA/dms"
cd "$PROJECT_DATA/dms"

# BRCA1 - Findlay 2018 SGE (MaveDB urn:mavedb:00000003)
echo "Downloading BRCA1 DMS (Findlay 2018)..."
mkdir -p brca1
wget -q -O brca1/findlay2018_scores.csv \
    "https://www.mavedb.org/api/v1/score-sets/urn:mavedb:00000003-a-1/scores/" 2>/dev/null || \
    echo "  NOTE: MaveDB API may need manual download. Check https://www.mavedb.org/#/experiment-sets/urn:mavedb:00000003"

# TP53 - Giacomelli 2018 (MaveDB urn:mavedb:00000068-b-1)
echo "Downloading TP53 DMS (Giacomelli 2018)..."
mkdir -p tp53
wget -q -O tp53/giacomelli2018_scores.csv \
    "https://www.mavedb.org/api/v1/score-sets/urn:mavedb:00000068-b-1/scores/" 2>/dev/null || \
    echo "  NOTE: MaveDB API may need manual download. Check https://www.mavedb.org/#/experiment-sets/urn:mavedb:00000068"

# CHEK2 - McCarthy-Leo 2024
echo "CHEK2 DMS (McCarthy-Leo 2024)..."
mkdir -p chek2
echo "  NOTE: Check MaveDB or paper supplementary for CHEK2 DMS data."
echo "  Paper: McCarthy-Leo 2024, 4,885 SNVs"

# DNMT3A - Huang 2022 (Blood) / Huang 2024
echo "DNMT3A DMS (Huang 2022/2024)..."
mkdir -p dnmt3a
echo "  NOTE: Check MaveDB for DNMT3A paired WT/R882H DMS data."
echo "  Paper: Huang 2022 (Blood), 4,020 single-AA substitutions"

# -------------------------------------------------------------------------
# 5. TERT Promoter MPRA (Kircher 2019)
# -------------------------------------------------------------------------
echo ""
echo "=== 5. TERT Promoter MPRA (project-specific) ==="
mkdir -p "$PROJECT_DATA/mpra"
cd "$PROJECT_DATA/mpra"

# MaveDB urn:mavedb:00000024-a or GEO GSE126550
echo "Downloading TERT promoter MPRA (Kircher 2019)..."
wget -q -O tert_mpra_scores.csv \
    "https://www.mavedb.org/api/v1/score-sets/urn:mavedb:00000024-a-1/scores/" 2>/dev/null || \
    echo "  NOTE: MaveDB API may need manual download. Check MaveDB urn:mavedb:00000024-a"
echo "  Also available from GEO: GSE126550"

# -------------------------------------------------------------------------
# 6. Benchmark Tool Precomputed Scores
# -------------------------------------------------------------------------
echo ""
echo "=== 6. Benchmark Tool Scores ==="
mkdir -p "$BASE/benchmarks"
cd "$BASE/benchmarks"

# CADD v1.7 (hg38 native) - whole genome precomputed
echo "Downloading CADD v1.7 (hg38)..."
mkdir -p cadd
cd cadd
if [ ! -f "whole_genome_SNVs.tsv.gz" ]; then
    echo "  Downloading whole-genome CADD scores (80 GB)..."
    echo "  NOTE: This is very large. Consider downloading only relevant chromosomes."
    echo "  URL: https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
    echo "  For now, we'll use the CADD web API for targeted scoring."
    # Uncomment below for full download:
    # wget -q "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
    # wget -q "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi"
fi
cd "$BASE/benchmarks"

# AlphaMissense precomputed
echo "Downloading AlphaMissense scores..."
mkdir -p alphamissense
cd alphamissense
if [ ! -f "AlphaMissense_hg38.tsv.gz" ]; then
    echo "  URL: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
    wget -q "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz" 2>/dev/null || \
        echo "  NOTE: May need manual download from Zenodo/Google Storage"
fi
cd "$BASE/benchmarks"

# SpliceAI precomputed
echo "SpliceAI precomputed scores..."
mkdir -p spliceai
echo "  NOTE: Download from Illumina BaseSpace or use SpliceAI Python package"
echo "  pip install spliceai; spliceai uses keras for real-time scoring"

# REVEL
echo "REVEL precomputed scores..."
mkdir -p revel
cd revel
if [ ! -f "revel_all_chromosomes.csv.zip" ]; then
    echo "  Downloading REVEL scores..."
    wget -q "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" -O revel_all_chromosomes.csv.zip 2>/dev/null || \
        echo "  NOTE: Download from https://sites.google.com/site/revelgenomics/downloads"
fi
cd "$BASE/benchmarks"

# PROVEAN (for indel comparison)
echo "PROVEAN..."
echo "  NOTE: Use PROVEAN web server or local installation"
echo "  http://provean.jcvi.org/index.php"

# -------------------------------------------------------------------------
# 7. Non-coding Benchmark Tools (hg19 — need liftover to hg38)
# -------------------------------------------------------------------------
echo ""
echo "=== 7. Non-coding Tools (hg19 → liftover needed) ==="

mkdir -p "$BASE/noncoding"
cd "$BASE/noncoding"

# LINSIGHT
echo "LINSIGHT scores (hg19)..."
mkdir -p linsight
echo "  URL: http://compgen.cshl.edu/LINSIGHT/LINSIGHT.bw"
echo "  Download the BigWig file and liftover to hg38"

# Eigen
echo "Eigen scores (hg19)..."
mkdir -p eigen
echo "  URL: https://xioniti01.u.hpc.mssm.edu/v1.1/"
echo "  Download per-chromosome Eigen scores; liftover via ANNOVAR"

# ncER
echo "ncER scores (hg19)..."
mkdir -p ncer
echo "  URL: https://github.com/TelentiLab/ncER_datasets"
echo "  Download percentile scores; liftover needed"

# Liftover chain file
echo "Downloading hg19 → hg38 liftover chain..."
mkdir -p "$BASE/liftover"
cd "$BASE/liftover"
if [ ! -f "hg19ToHg38.over.chain.gz" ]; then
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
fi

# -------------------------------------------------------------------------
# 8. ENCODE ChIP-seq Peaks (hg38)
# -------------------------------------------------------------------------
echo ""
echo "=== 8. ENCODE Regulatory Data ==="
mkdir -p "$BASE/encode"
cd "$BASE/encode"

echo "ENCODE data will be downloaded per-gene region from ENCODE portal."
echo "Key datasets: DHS, H3K27ac, H3K4me1, CTCF ChIP-seq"
echo "Use ENCODE REST API: https://www.encodeproject.org/help/rest-api/"

# -------------------------------------------------------------------------
# 9. gnomAD Constraint Data
# -------------------------------------------------------------------------
echo ""
echo "=== 9. gnomAD Constraint ==="
mkdir -p "$BASE/gnomad"
cd "$BASE/gnomad"

echo "Downloading gnomAD constraint metrics..."
# Gene-level constraint (LOEUF)
if [ ! -f "gnomad.v4.1.constraint_metrics.tsv" ]; then
    wget -q "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" 2>/dev/null || \
        echo "  NOTE: Check gnomAD downloads page for latest constraint file"
fi
echo "  Also need: Gnocchi regional constraint (1kb resolution)"
echo "  URL: https://gnomad.broadinstitute.org/downloads#v4-constraint"

# -------------------------------------------------------------------------
# 10. GRCm39 Mouse Genome (for cross-species analysis)
# -------------------------------------------------------------------------
echo ""
echo "=== 10. GRCm39 Mouse Genome ==="
mkdir -p "$BASE/mouse"
cd "$BASE/mouse"

if [ ! -f "GRCm39.primary_assembly.genome.fa" ]; then
    echo "Downloading GRCm39 mouse genome..."
    wget -q "https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    mv Mus_musculus.GRCm39.dna.primary_assembly.fa GRCm39.primary_assembly.genome.fa
    echo "Building index..."
    python3 -c "from pyfaidx import Fasta; f = Fasta('GRCm39.primary_assembly.genome.fa'); print(f'Indexed: {len(f.keys())} contigs')"
else
    echo "Mouse genome already exists, skipping."
fi

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
echo ""
echo "=== Download Summary ==="
echo "Data directory: $BASE"
du -sh "$BASE"/* 2>/dev/null || true
echo ""
echo "Manual downloads still needed:"
echo "  1. DMS data: verify MaveDB API downloads or get from paper supplements"
echo "  2. CADD v1.7: full genome (~80GB) or use web API for targeted scoring"
echo "  3. LINSIGHT, Eigen, ncER: download + liftover hg19→hg38"
echo "  4. ENCODE: per-gene ChIP-seq peaks via REST API"
echo "  5. gnomAD Gnocchi: regional constraint at 1kb resolution"
echo ""
echo "=== Done: $(date) ==="
