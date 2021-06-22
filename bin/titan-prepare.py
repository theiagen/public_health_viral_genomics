#! /usr/bin/env python3
"""
usage: titan-prepare [-h] [-f STR] [--fastq_seperator STR] [--fastq_pattern STR] [--pe1_pattern STR] [--pe2_pattern STR] 
                     [-r] [--prefix STR] FASTQ_PATH RUN_ID PLATFORM

titan-prepare - Read a directory and prepare a JSON for input to Titan

positional arguments:
  FASTQ_PATH            Directory where FASTQ files are stored
  RUN_ID                Run ID to associate with the samples.
  PLATFORM              The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont

optional arguments:
  -h, --help            show this help message and exit
  -f STR, --fastq_ext STR
                        Extension of the FASTQs. Default: .fastq.gz
  --fastq_seperator STR
                        Split FASTQ name on the last occurrence of the separator. Default: _
  --fastq_pattern STR   Glob pattern to match FASTQs. Default: *.fastq.gz
  --pe1_pattern STR     Designates difference first set of paired-end reads. Default: ([Aa]|[Rr]1|1) (R1, r1, 1, A, a)
  --pe2_pattern STR     Designates difference second set of paired-end reads. Default: ([Bb]|[Rr]2|2) (R2, r2, 2, AB b)
  -r, --recursive       Directories will be traversed recursively
  --prefix STR          Replace the absolute path with a given string. Default: Use absolute path

Heavily based off bactopia prepare (https://github.com/bactopia/bactopia/blob/master/bin/helpers/bactopia-prepare.py)
"""

def search_path(path, pattern, recursive=False):
    from pathlib import Path
    if recursive:
        return Path(path).rglob(pattern)
    else:
        return Path(path).glob(pattern)


def get_path(fastq, abspath, prefix):
    fastq_path = str(fastq.absolute())
    if prefix:
        return fastq_path.replace(abspath, prefix.rstrip("/"))
    return fastq_path


if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import glob
    import json
    import os
    import re
    import sys

    parser = ap.ArgumentParser(
        prog='titan-prepare',
        conflict_handler='resolve',
        description=(
            f'titan-prepare - Read a directory and prepare a JSON for input to Titan'
        )
    )
    parser.add_argument('path', metavar="FASTQ_PATH", type=str,
                        help='Directory where FASTQ files are stored')
    parser.add_argument('run_id', metavar='RUN_ID', type=str, help='Run ID to associate with the samples.')
    parser.add_argument(
        'platform', metavar='PLATFORM', type=str, choices=['clearlabs', 'illumina_pe', 'illumina_se', 'ont'],
        help='The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont'
    )
    parser.add_argument(
        '--primers', metavar='STR', type=str, default="",
        help='A file containing primers (bed format) used during sequencing.'
    )
    parser.add_argument(
        '-f', '--fastq_ext', metavar='STR', type=str, default=".fastq.gz",
        help='Extension of the FASTQs. Default: .fastq.gz'
    )
    parser.add_argument(
        '--fastq_seperator', metavar='STR', type=str, default="_",
        help='Split FASTQ name on the last occurrence of the separator. Default: _'
    )
    parser.add_argument(
        '--fastq_pattern', metavar='STR', type=str, default="*.fastq.gz",
        help='Glob pattern to match FASTQs. Default: *.fastq.gz'
    )
    parser.add_argument(
        '--pe1_pattern', metavar='STR', type=str, default="[Aa]|[Rr]1|1",
        help='Designates difference first set of paired-end reads. Default: ([Aa]|[Rr]1|1) (R1, r1, 1, A, a)'
    )
    parser.add_argument(
        '--pe2_pattern', metavar='STR', type=str, default="[Bb]|[Rr]2|2",
        help='Designates difference second set of paired-end reads. Default: ([Bb]|[Rr]2|2) (R2, r2, 2, AB b)'
    )
    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Directories will be traversed recursively'
    )
    parser.add_argument(
        '--prefix', metavar='STR', type=str,
        help='Replace the absolute path with a given string. Default: Use absolute path'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    abspath = os.path.abspath(args.path)
    SAMPLES = {}

    # Match FASTQS
    for fastq in search_path(abspath, args.fastq_pattern, recursive=args.recursive):
        fastq_name = fastq.name.replace(args.fastq_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(args.fastq_seperator, 1)
        sample_name = split_vals[0]
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {'pe': {'r1': [], 'r2': []}, 'se': []}

        if len(split_vals) == 1:
            # single-end
            SAMPLES[sample_name]['se'].append(get_path(fastq, abspath, args.prefix))
        else:
            # paired-end
            pe1 = re.compile(args.pe1_pattern)
            pe2 = re.compile(args.pe2_pattern)
            if pe1.match(split_vals[1]):
                SAMPLES[sample_name]['pe']['r1'].append(get_path(fastq, abspath, args.prefix))
            elif pe2.match(split_vals[1]):
                SAMPLES[sample_name]['pe']['r2'].append(get_path(fastq, abspath, args.prefix))
            else:
                print(f'ERROR: Could not determine read set for "{fastq_name}".', file=sys.stderr)
                print(f'ERROR: Found {split_vals[1]} expected (R1: {args.pe1_pattern} or R2: {args.pe2_pattern})', file=sys.stderr)
                print(f'ERROR: Please use --pe1_pattern and --pe2_pattern to correct and try again.', file=sys.stderr)
                sys.exit(1)

    FOFN = []
    for sample, vals in sorted(SAMPLES.items()):
        r1_reads = vals['pe']['r1']
        r2_reads = vals['pe']['r2']
        se_reads = vals['se']
        errors = []
        pe_count = len(r1_reads) + len(r2_reads)

        # Validate everything
        if len(r1_reads) != len(r2_reads):
            # PE reads must be a pair
            errors.append(f'ERROR: "{sample}" must have equal paired-end read sets (R1 has {len(r1_reads)} and R2 has {len(r2_reads)}, please check.')
        elif pe_count > 2:
            # PE reads must be a pair
            errors.append(f'ERROR: "{sample}" cannot have more than two paired-end FASTQ, please check.')

        if len(se_reads) > 1:
            # Can't have multiple SE reads
            errors.append(f'ERROR: "{sample}" has more than two single-end FASTQs, please check.')
        elif pe_count and len(se_reads):
            # Can't have SE and PE reads unless long reads
            errors.append(f'ERROR: "{sample}" has paired and single-end FASTQs, please check.')

        if errors:
            print('\n'.join(errors), file=sys.stderr)
        else:
            r1 = ''
            r2 = ''

            if pe_count:
                platform = 'illumina_pe'
                r1 = r1_reads[0]
                r2 = r2_reads[0]

            if se_reads:
                r1 = se_reads[0]

            FOFN.append({
                'samplename': sample, 
                'run_id': args.run_id,
                'platform': args.platform,
                'r1': r1,
                'r2': r2
            })

    if FOFN:
        inputs_json = {
            "titan_gc.samples": FOFN,
            "titan_gc.primer_bed": args.primers
        }
        print(json.dumps(inputs_json, indent = 4))
