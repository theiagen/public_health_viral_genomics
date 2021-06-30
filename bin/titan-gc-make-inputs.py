#! /usr/bin/env python3
"""
usage: titan-gc-make-inputs [-h] [--sample STR] [--run_id STR] [--platform STR] [--r1 STR] [--r2 STR] [--primers PRIMER]

titan-gc-make-inputs - Creates a input JSON for a single sample

optional arguments:
  -h, --help        show this help message and exit
  --sample STR      Directory where FASTQ files are stored
  --run_id STR      Run ID to associate with the samples.
  --platform STR    The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont
  --r1 STR          R1 FASTQ of a read pair or single-end FASTQ.
  --r2 STR          R2 FASTQ of a read pair
  --primers PRIMER  A file containing primers (bed format) used during sequencing.
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

    parser.add_argument(
        '--sample', metavar="STR", type=str, help='Directory where FASTQ files are stored'
    )
    parser.add_argument(
        '--run_id', metavar='STR', type=str, help='Run ID to associate with the samples.'
    )
    parser.add_argument(
        '--platform', metavar='STR', type=str, choices=['clearlabs', 'illumina_pe', 'illumina_se', 'ont'],
        help='The platform used for sequencing. Options: clearlabs, illumina_pe, illumina_se, ont'
    )
    parser.add_argument(
        '--r1', metavar='STR', type=str,
        help='R1 FASTQ of a read pair or single-end FASTQ.'
    )
    parser.add_argument(
        '--r2', metavar='STR', type=str,
        help='R2 FASTQ of a read pair'
    )
    parser.add_argument(
        '--primers', metavar='PRIMER', type=str,
        help='A file containing primers (bed format) used during sequencing.'
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
    print(json.dumps(inputs, indent = 4))
