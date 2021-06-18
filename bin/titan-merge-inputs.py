#! /usr/bin/env python3
"""
usage: titan-merge-inputs [-h] [inputs [inputs ...]]

titan-merge-inputs - Merge multiple JSON files produced Titan prepare

positional arguments:
  inputs

optional arguments:
  -h, --help  show this help message and exit
"""

if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import json
    import sys

    parser = ap.ArgumentParser(
        prog='titan-merge-inputs',
        conflict_handler='resolve',
        description=(
            f'titan-merge-inputs - Merge multiple JSON files produced Titan prepare'
        )
    )
    parser.add_argument('inputs', nargs='*')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    json_data = []
    for j in args.inputs:
        with open(j, 'rt') as json_fh:
            json_data.append(json.load(json_fh))

    primers = []
    samples = []
    for j in json_data:
        if j["titan.primer_bed"]:
            primers.append(j["titan.primer_bed"])
        samples.extend(j["titan.samples"])

    primers = list(set(primers))
    primer_bed = ""
    if len(primers) == 1:
        primer_bed = primers[0]
    elif len(primers) >1:
        print(f'Multiple primer files given, can only have one.', file=sys.stderr)
        for p in primers:
            print(f'\tFound: {p}', file=sys.stderr)
        sys.exit(1)
    
    merged_inputs = {
        "titan.samples": samples,
        "titan.primer_bed": primer_bed
    }
    print(json.dumps(merged_inputs, indent=4))
