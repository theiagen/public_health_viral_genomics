#! /usr/bin/env python3
"""
usage: titan-prepare [-h] [--outdir STR] [--debug] METADATA_JSON

titan-prepare - Read a directory and prepare a JSON for input to Titan

positional arguments:
  METADATA_JSON  The metadata.json output (-m) from the Titan run.

optional arguments:
  -h, --help     show this help message and exit
  --outdir STR   Directory to copy files to. (Default: ./titan).
  --debug        Print helpful information
"""

def mkdir(path):
    from pathlib import Path
    Path(path).mkdir(parents=True, exist_ok=True)

def read_titan_results(tsv):
    results = {}
    with open(tsv, 'rt') as tsv_fh:
        for line in tsv_fh:
            if "result-header" not in results:
                results["result-header"] = line
            else:
                # store the line, col[0] is samplename
                results[line.split("\t")[0]] = line
    return results

if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    from shutil import copy2
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
    parser.add_argument('metadata', metavar="METADATA_JSON", type=str,
                        help='The metadata.json output (-m) from the Titan run.')
    parser.add_argument('--outdir', metavar='STR', type=str, default="./titan",
                        help='Directory to copy files to. (Default: ./titan).')
    parser.add_argument('--debug', action='store_true', help='Print helpful information')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if not os.path.exists(args.metadata):
        print(f"Unable to find {args.metadata}, please verify the path.", file=sys.stderr)
        sys.exit(1)
    
    metadata = None
    with open(args.metadata, 'rt') as metadata_fh:
        metadata = json.load(metadata_fh)
    """
    Read sample inputs

    The samples included in the run are found at metadata["inputs"]["samples"] (list of dicts)

    Example print(metadata["inputs"]["samples"])
    [{
      "samplename": "sample_01",
      "run_id": "run_01",
      "r1": "sample_01.fastq.gz",
      "r2": "",
      "platform": "clearlabs"
    }, {
      "samplename": "sample_02",
      "run_id": "run_02",
      "r1": "sample_02_R1.fastq.gz",
      "r2": "sample_02_R2.fastq.gz",
      "platform": "illumina_pe"
    }]
    """
    samples = {}
    runs = {}
    for sample in metadata["inputs"]["samples"]:
        samples[sample["samplename"]] = sample
        if sample["run_id"] not in runs:
            runs[sample["run_id"]] = []
        runs[sample["run_id"]].append(sample["samplename"])

    """
    Start moving files: metadata["outputs"]
    Output keys:
        "titan.reads_dehosted"
        "titan.kraken_report"
        "titan.pango_lineage_report"
        "titan.nextclade_json"
        "titan.consensus_flagstat"
        "titan.kraken_report_dehosted"
        "titan.vadr_alerts_list"
        "titan.aligned_bam"
        "titan.assembly_fasta"
        "titan.nextclade_tsv"
        "titan.consensus_stats"
        "titan.auspice_json"
        "titan.aligned_bai"

    The merged summary is under: "titan.merged_summaries"
    Output files should start with the sample name (sample01.file or sample01_file)
    """
    mkdir(f"{args.outdir}")
    titan_results = None
    for key, outputs in metadata["outputs"].items():
        task_name = key.replace("titan.", "")
        if args.debug:
            print(f"Working on {task_name} outputs", file=sys.stderr)
        if key == "titan.merged_summaries":
            print(outputs)
            if args.debug:
                print(f"Copying {outputs} to {args.outdir}/complete-titan-results.tsv", file=sys.stderr)
            copy2(outputs, f"{args.outdir}/titan-results.tsv")
            titan_results = read_titan_results(outputs)
        else:
            for output in outputs:
                filename = os.path.basename(output)
                samplename = filename.split(".")[0]
                if samplename not in samples:
                    samplename = filename.split("_")[0]
                if samplename not in samples:
                    print(f"Unable to associate {output} with a sample", file=sys.stderr)
                    continue
                
                copy_path = f"{args.outdir}/{samples[samplename]['run_id']}/{task_name}"
                mkdir(copy_path)

                if args.debug:
                    print(f"Copying {output} to {copy_path}", file=sys.stderr)
                copy2(output, copy_path)

    # Write per-run Titan results
    for run_id, inputs in runs.items():
        run_results = f"{args.outdir}/{run_id}/{run_id}-titan-results.tsv"
        if args.debug:
            print(f"Writing Titan results for {run_id} to {run_results}", file=sys.stderr)
        with open(run_results, 'wt') as result_fh:
            result_fh.write(titan_results["result-header"])
            for sample in inputs:
                 result_fh.write(titan_results[sample])
