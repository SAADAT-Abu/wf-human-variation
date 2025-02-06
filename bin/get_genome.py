#!/usr/bin/env python
"""Check BAM header and return genome build based on chromosome sizes."""

import argparse
import os
import sys


HG38_URL = (
    "https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/ref/" +
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
CHROMOSOME_SIZES = {
    'hg19': {
        'chr1': '249250621',
        'chr2': '243199373',
        'chr3': '198022430',
        'chr4': '191154276',
        'chr5': '180915260',
        'chr6': '171115067',
        'chr7': '159138663',
        'chr8': '146364022',
        'chr9': '141213431',
        'chr10': '135534747',
        'chr11': '135006516',
        'chr12': '133851895',
        'chr13': '115169878',
        'chr14': '107349540',
        'chr15': '102531392',
        'chr16': '90354753',
        'chr17': '81195210',
        'chr18': '78077248',
        'chr19': '59128983',
        'chr20': '63025520',
        'chr21': '48129895',
        'chr22': '51304566',
        'chrX': '155270560',
        'chrY': '59373566'},
    'hg38': {
        'chr1': '248956422',
        'chr2': '242193529',
        'chr3': '198295559',
        'chr4': '190214555',
        'chr5': '181538259',
        'chr6': '170805979',
        'chr7': '159345973',
        'chr8': '145138636',
        'chr9': '138394717',
        'chr10': '133797422',
        'chr11': '135086622',
        'chr12': '133275309',
        'chr13': '114364328',
        'chr14': '107043718',
        'chr15': '101991189',
        'chr16': '90338345',
        'chr17': '83257441',
        'chr18': '80373285',
        'chr19': '58617616',
        'chr20': '64444167',
        'chr21': '46709983',
        'chr22': '50818468',
        'chrX': '156040895',
        'chrY': '57227415'},
    'T2T': {
        'NC_060925.1': '248387328',
        'NC_060926.1': '242696752',
        'NC_060927.1': '201105948',
        'NC_060928.1': '193574945',
        'NC_060929.1': '182045439',
        'NC_060930.1': '172126628',
        'NC_060931.1': '160567428',
        'NC_060932.1': '146259331',
        'NC_060933.1': '150617247',
        'NC_060934.1': '134758134',
        'NC_060935.1': '135127769',
        'NC_060936.1': '133324548',
        'NC_060937.1': '113566686',
        'NC_060938.1': '101161492',
        'NC_060939.1': '99753195',
        'NC_060940.1': '96330374',
        'NC_060941.1': '84276897',
        'NC_060942.1': '80542538',
        'NC_060943.1': '61707364',
        'NC_060944.1': '66210255',
        'NC_060945.1': '45090682',
        'NC_060946.1': '51324926',
        'NC_060947.1': '154259566',
        'NC_060948.1': '62460029'}
}



ALLOWED_CHR = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", 
    "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1",
    "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1",
    "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1",
    "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1",
    "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"
]


def chromosome_sizes(extracted_sizes):
    """Get dictionary of chromosomes and sizes from samtools index."""
    fasta_sizes = {}
    with open(extracted_sizes, 'r') as fa_idx:
        for line in fa_idx:
            line = line.rstrip()
            cols = line.split('\t')
            # Do not prepend 'chr' for T2T chromosomes
            if cols[0] in CHROMOSOME_SIZES['T2T'] or cols[0].startswith("NC_"):
                fasta_sizes[cols[0]] = cols[1]  # Ensure string format
            elif cols[0] in ALLOWED_CHR:
                fasta_sizes[cols[0]] = cols[1]
    return fasta_sizes


def get_genome(sizes):
    """Get genome based on chromosome sizes."""
    if not sizes:
        return ""
    for known_genome_build in CHROMOSOME_SIZES.keys():
        if sizes == CHROMOSOME_SIZES[known_genome_build]:  # Ensure exact match
            return known_genome_build
    return ""

def check_genome(genome_build, str_flag):
    """Determine if genome is suitable for this workflow."""
    bad_genome = False
    extra_msg_context = (
        "#####################################################################\n"
        "# INPUT DATA PROBLEM\n"
        "The genome build detected in the BAM is not compatible with this\n"
        "workflow as it does not appear to be hg19/GRCh37, hg38/GRCh38 or T2TCHM13v2.\n"
        "If you are trying to run this workflow with non-human data, please\n"
        "consult the 'Genome compatibility and running the workflow on\n"
        "non-human genomes' section of the README.\n"
        "####################################################################\n"
    )
    if not genome_build:
        bad_genome = True
    elif str_flag and genome_build != "hg38":
        bad_genome = True
        extra_msg_context = (
            "#####################################################################\n"
            "# INPUT DATA PROBLEM\n"
            "The genome build detected in the BAM is not compatible with this\n"
            "workflow.\n"
            f"Detected genome: {genome_build}, but genotyping STRs can only be\n"
            "performed when aligned to build 38.\n"
            "To perform STR calling, you need to run the workflow providing the\n"
            "following reference genome to the --ref parameter:\n\n"
            f"{HG38_URL}\n\n"
            "Alternatively, disable STR calling by setting --str false.\n"
            "####################################################################\n"
        )
    return (bad_genome, extra_msg_context)


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--chr_counts', required=True, dest="chr_counts",
        help="Output from samtools faidx"
    )
    parser.add_argument(
        '-o', '--output', required=True, dest="output",
        help="Output genome"
    )
    parser.add_argument(
        '--str',  action='store_true',
        dest="str_",
        default=False,
        help="STR flag"
    )
    args = parser.parse_args()

    all_sizes = chromosome_sizes(args.chr_counts)

    genome_build = get_genome(all_sizes)
    bad_genome, extra_msg_context = check_genome(genome_build, args.str_)

    # explode on bad genome
    if bad_genome:
        sys.stderr.write(extra_msg_context)
        sys.exit(os.EX_DATAERR)

    # otherwise write out the genome name
    result = open(args.output, 'w')
    result.write(genome_build)


if __name__ == '__main__':
    main()
