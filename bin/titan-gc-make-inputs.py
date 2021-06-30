#! /usr/bin/env python3
"""
usage: titan-gc-make-inputs [-h] [--sample STR] [--run_id STR] [--platform STR] [--r1 STR] [--r2 STR] [--primers PRIMER] 
                            [--pangolin_docker STR] [--clearlabs_normalise INT] [--ont_normalise INT] [--seq_method STR]

titan-gc-make-inputs - Creates a input JSON for a single sample

optional arguments:
  -h, --help            show this help message and exit

Titan-GC Prepare Parameters:
  --sample STR          Directory where FASTQ files are stored
  --run_id STR          Run ID to associate with the samples.
  --platform STR        The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont
  --r1 STR              R1 FASTQ of a read pair or single-end FASTQ.
  --r2 STR              R2 FASTQ of a read pair
  --primers PRIMER      A file containing primers (bed format) used during sequencing.

Optional Titan-GC Workflow Parameters:
  --pangolin_docker STR
                        Docker image used to run Pangolin
  --clearlabs_normalise INT
                        Value to normalize Clearlabs read counts
  --ont_normalise INT   Value to normalize ONT read counts
  --seq_method STR      Seqeuncing method used
"""

if __name__ == '__main__':
    import argparse as ap
    from pathlib import Path
    import json
    import os
    import sys

    parser = ap.ArgumentParser(
        prog='titan-gc-make-inputs',
        conflict_handler='resolve',
        description=(
            f'titan-gc-make-inputs - Creates a input JSON for a single sample'
        )
    )

    group1 = parser.add_argument_group('Titan-GC Prepare Parameters')
    group1.add_argument(
        '--sample', metavar="STR", type=str, help='Directory where FASTQ files are stored'
    )
    group1.add_argument(
        '--run_id', metavar='STR', type=str, help='Run ID to associate with the samples.'
    )
    group1.add_argument(
        '--platform', metavar='STR', type=str, choices=['clearlabs', 'illumina_pe', 'illumina_se', 'ont'],
        help='The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont'
    )
    group1.add_argument(
        '--r1', metavar='STR', type=str,
        help='R1 FASTQ of a read pair or single-end FASTQ.'
    )
    group1.add_argument(
        '--r2', metavar='STR', type=str,
        help='R2 FASTQ of a read pair'
    )
    group1.add_argument(
        '--primers', metavar='PRIMER', type=str,
        help='A file containing primers (bed format) used during sequencing.'
    )

    group2 = parser.add_argument_group('Optional Titan-GC Workflow Parameters')
    group2.add_argument(
        '--pangolin_docker', metavar='STR', type=str,
        help='Docker image used to run Pangolin'
    )
    group2.add_argument(
        '--clearlabs_normalise', metavar='INT', type=str,
        help='Value to normalize Clearlabs read counts'
    )
    group2.add_argument(
        '--ont_normalise', metavar='INT', type=str,
        help='Value to normalize ONT read counts'
    )
    group2.add_argument(
        '--seq_method', metavar='STR', type=str, help='Seqeuncing method used'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    EMPTY_FASTQ = f"{str(Path.home())}/.titan/EMPTY.fastq.gz"
    if not os.path.exists(EMPTY_FASTQ):
        Path(f"{str(Path.home())}/.titan").mkdir(parents=True, exist_ok=True)
        with open(EMPTY_FASTQ, 'a'):
            pass

    inputs = {
        "titan_gc.samples": [
            {
                'samplename': args.sample, 
                'run_id': args.run_id,
                'platform': args.platform,
                'r1': str(Path(args.r1).absolute()),
                'r2': str(Path(args.r2).absolute()) if args.r2 else EMPTY_FASTQ,
                'primer_bed': str(Path(args.primers).absolute())
            }
        ]
    }

    # Add optional parameters if user specified them
    if args.pangolin_docker:
        inputs['titan_gc.pangolin_docker_image'] = args.pangolin_docker
    if args.clearlabs_normalise:
        inputs['titan_gc.titan_clearlabs.normalise'] = args.clearlabs_normalise
    if args.ont_normalise:
        inputs['titan_gc.titan_ont.normalise'] = args.ont_normalise
    if args.seq_method:
        key = {'clearlabs': "titan_gc.titan_clearlabs.seq_method", 'illumina_pe': "titan_gc.titan_illumina_pe.seq_method",
                'illumina_se': "titan_gc.titan_illumina_se.seq_method", 'ont': "titan_gc.titan_ont.seq_method"}
        inputs[key[args.platform]] = args.seq_method
    
    print(json.dumps(inputs, indent = 4))
