#! /usr/bin/env python3
"""
usage: titan-gc-prepare [-h] [-f STR] [--fastq_separator STR] [--fastq_pattern STR] [--pe1_pattern STR] [--pe2_pattern STR]
                        [-r] [--prefix STR] [--tsv] [--pangolin_docker STR] [--clearlabs_normalise INT] [--ont_normalise INT]
                        [--seq_method STR] FASTQ_PATH WORKFLOW PRIMER

titan-gc-prepare - Read a directory and prepare a JSON for input to Titan GC

optional arguments:
  -h, --help            show this help message and exit

Titan-GC Prepare Parameters:
  FASTQ_PATH            Directory where FASTQ files are stored
  WORKFLOW              The TItan-GC workflow to use for anlaysis. Options: clearlabs, illumina_pe, illumina_se, ont
  PRIMERS               A file containing primers (bed format) used during sequencing.
  -f STR, --fastq_ext STR
                        Extension of the FASTQs. Default: .fastq.gz
  --fastq_separator STR
                        Split FASTQ name on the last occurrence of the separator. Default: _
  --fastq_pattern STR   Glob pattern to match FASTQs. Default: *.fastq.gz
  --pe1_pattern STR     Designates difference first set of paired-end reads. Default: ([Aa]|[Rr]1|1) (R1, r1, 1, A, a)
  --pe2_pattern STR     Designates difference second set of paired-end reads. Default: ([Bb]|[Rr]2|2) (R2, r2, 2, AB b)
  -r, --recursive       Directories will be traversed recursively
  --prefix STR          Replace the absolute path with a given string. Default: Use absolute path
  --tsv                 Output FOFN as a TSV (Default JSON)

Optional Titan-GC Workflow Parameters:
  --pangolin_docker STR
                        Docker image used to run Pangolin
  --clearlabs_normalise INT
                        Value to normalize Clearlabs read counts
  --ont_normalise INT   Value to normalize ONT read counts
  --seq_method STR      Seqeuncing method used

Heavily based off bactopia prepare (https://github.com/bactopia/bactopia/blob/master/bin/helpers/bactopia-prepare.py)
"""
from pathlib import Path

def search_path(path, pattern, recursive=False):
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
        prog='titan-gc-prepare',
        conflict_handler='resolve',
        description=(
            f'titan-gc-prepare - Read a directory and prepare a JSON for input to Titan GC'
        )
    )
    group1 = parser.add_argument_group('Titan-GC Prepare Parameters')
    group1.add_argument('path', metavar="FASTQ_PATH", type=str,
                        help='Directory where FASTQ files are stored')
    group1.add_argument(
        'workflow', metavar='WORKFLOW', type=str, choices=['clearlabs', 'illumina_pe', 'illumina_se', 'ont'],
        help='The Titan-GC workflow to use for analysis. Options: clearlabs, illumina_pe, illumina_se, ont'
    )
    group1.add_argument(
        'primers', metavar='PRIMER', type=str, default="",
        help='A file containing primers (bed format) used during sequencing.'
    )
    group1.add_argument(
        '-f', '--fastq_ext', metavar='STR', type=str, default=".fastq.gz",
        help='Extension of the FASTQs. Default: .fastq.gz'
    )
    group1.add_argument(
        '--fastq_separator', metavar='STR', type=str, default="_",
        help='Split FASTQ name on the last occurrence of the separator. Default: _'
    )
    group1.add_argument(
        '--fastq_pattern', metavar='STR', type=str, default="*.fastq.gz",
        help='Glob pattern to match FASTQs. Default: *.fastq.gz'
    )
    group1.add_argument(
        '--pe1_pattern', metavar='STR', type=str, default="[Aa]|[Rr]1|1",
        help='Designates difference first set of paired-end reads. Default: ([Aa]|[Rr]1|1) (R1, r1, 1, A, a)'
    )
    group1.add_argument(
        '--pe2_pattern', metavar='STR', type=str, default="[Bb]|[Rr]2|2",
        help='Designates difference second set of paired-end reads. Default: ([Bb]|[Rr]2|2) (R2, r2, 2, AB b)'
    )
    group1.add_argument(
        '-r', '--recursive', action='store_true',
        help='Directories will be traversed recursively'
    )
    group1.add_argument(
        '--prefix', metavar='STR', type=str,
        help='Replace the absolute path with a given string. Default: Use absolute path'
    )
    group1.add_argument('--tsv', action='store_true', help='Output FOFN as a TSV (Default JSON)')

    group2 = parser.add_argument_group('Optional Titan-GC Workflow Parameters')
    group2.add_argument(
        '--pangolin_docker', metavar='STR', type=str,
        help='Docker image used to run Pangolin (takes priority over --params)'
    )
    group2.add_argument(
        '--params', metavar='STR', type=str, help='A JSON file containing parameter values'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    abspath = os.path.abspath(args.path)
    SAMPLES = {}
    EMPTY_FASTQ = f"{str(Path.home())}/.titan/EMPTY.fastq.gz"
    if not os.path.exists(EMPTY_FASTQ):
        Path(f"{str(Path.home())}/.titan").mkdir(parents=True, exist_ok=True)
        with open(EMPTY_FASTQ, 'a'):
            pass

    # Match FASTQS
    for fastq in search_path(abspath, args.fastq_pattern, recursive=args.recursive):
        fastq_name = fastq.name.replace(args.fastq_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(args.fastq_separator, 1)
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
                r1 = r1_reads[0]
                r2 = r2_reads[0]

            if se_reads:
                r1 = se_reads[0]
                r2 = EMPTY_FASTQ

            FOFN.append({
                'sample': sample, 
                'titan_wf': args.workflow,
                'r1': r1,
                'r2': r2,
                'primers': get_path(Path(args.primers), abspath, args.prefix)
            })

    if FOFN:
        if args.tsv:
            needs_header = True
            for f in FOFN:
                if needs_header:
                    print("\t".join(['sample', 'titan_wf', 'r1', 'r2', 'primers']))
                    needs_header = False
                print("\t".join([f['sample'], f['titan_wf'], f['r1'], f['r2'], f['primers']]))
        else:
            inputs_json = {
                "titan_gc.samples": FOFN
            }
            params_json = {}


            # Add optional parameters if user specified them
            if args.params:
                with open(args.params, 'rt') as json_fh:
                    params_json = json.load(json_fh)

            if args.pangolin_docker:
                params_json['titan_gc.pangolin_docker_image'] = args.pangolin_docker
            
            print(json.dumps({**inputs_json, **params_json}, indent = 4))
