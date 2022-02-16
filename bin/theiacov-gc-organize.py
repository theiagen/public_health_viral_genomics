#! /usr/bin/env python3
"""
usage: theiacov-gc-organize [-h] [--outdir STR] [--debug] METADATA_JSON

theiacov-gc-organize- Read Cromwell metadata to organize files into readable structure

positional arguments:
  METADATA_JSON  The metadata.json output (-m) from the TheiaCoV GC run.

optional arguments:
  -h, --help     show this help message and exit
  --outdir STR   Directory to copy files to. (Default: ./theiacov-gc).
  --debug        Print helpful information
"""
# Known file extensions from TheiaCoV-GC outputs
OUTPUTS = {
    'aligned_bam': {'ext': ['.primertrim.sorted.bam', '.primertrimmed.rg.sorted.bam'], 'folder': 'alignments'},
    'aligned_bai': {'ext': ['.primertrim.sorted.bam.bai', '.primertrimmed.rg.sorted.bam.bai'], 'folder': 'alignments'},
    'assembly_fasta': {'ext': ['.ivar.consensus.fasta', '.medaka.consensus.fasta'], 'folder': 'assemblies'},
    'auspice_json': {'ext': ['.ivar.consensus.nextclade.auspice.json', '.medaka.consensus.nextclade.auspice.json'], 'folder': 'nextclade'},
    'consensus_flagstat': {'ext': ['.flagstat.txt'], 'folder': 'alignments'},
    'consensus_stats': {'ext': ['.stats.txt'], 'folder': 'alignments'},
    'json_summary': {'ext': ['.results.json'], 'folder': 'results'},
    'nextclade_json': {'ext': ['.ivar.consensus.nextclade.json', '.medaka.consensus.nextclade.json'], 'folder': 'nextclade'},
    'nextclade_tsv': {'ext': ['.ivar.consensus.nextclade.tsv', '.medaka.consensus.nextclade.tsv'], 'folder': 'nextclade'},
    'pango_lineage_report': {'ext': ['.pangolin_report.csv'], 'folder': 'pangolin_reports'},
    'vadr_alerts_list': {'ext': ['.ivar.consensus.vadr.alt.list', '.medaka.consensus.vadr.alt.list'], 'folder': 'vadr_alerts'},
    'reads_dehosted': {'ext': ['_dehosted.fastq.gz'], 'folder': 'dehosted_reads'},
    'kraken_report_dehosted': {'ext': ['_kraken2_report.txt'], 'folder': 'kraken2_reports'},
    'kraken_report': {'ext': ['_kraken2_report.txt'], 'folder': 'kraken2_reports'},
    'consensus_variants': {'ext': ['.pass.vcf', '.variants.vcf'], 'folder': 'alignments'}
}
import json

def mkdir(path):
    from pathlib import Path
    Path(path).mkdir(parents=True, exist_ok=True)

def read_theiacov_results(tsv, is_json=False):
    results = {}
    with open(tsv, 'rt') as tsv_fh:
        for line in tsv_fh:
            if is_json:
                json_data = json.loads(line.rstrip())
                results[json_data['sample']] = line
            else:
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
    import os
    import re
    import sys

    parser = ap.ArgumentParser(
        prog='theiacov-gc-organize',
        conflict_handler='resolve',
        description=(
            f'theiacov-gc-organize- Read Cromwell metadata to organize files into readable structure'
        )
    )
    parser.add_argument('metadata', metavar="METADATA_JSON", type=str,
                        help='The metadata.json output (-m) from the TheiaCoV GC run.')
    parser.add_argument('--outdir', metavar='STR', type=str, default="./theiacov-gc",
                        help='Directory to copy files to. (Default: ./theiacov-gc).')
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
      "sample": "sample_01",
      "workflow": "clearlabs",
      "r1": "sample_01.fastq.gz",
      "r2": "EMPTY.fastq",
      "primers": "primers.bed"
    }]
    """
    samples = {}
    for sample in metadata["inputs"]["samples"]:
        samples[sample["sample"]] = sample

    """
    Start moving files: metadata["outputs"]
    Output keys:
        "theiacov_gc.reads_dehosted"
        "theiacov_gc.kraken_report"
        "theiacov_gc.pango_lineage_report"
        "theiacov_gc.nextclade_json"
        "theiacov_gc.consensus_flagstat"
        "theiacov_gc.kraken_report_dehosted"
        "theiacov_gc.vadr_alerts_list"
        "theiacov_gc.aligned_bam"
        "theiacov_gc.assembly_fasta"
        "theiacov_gc.nextclade_tsv"
        "theiacov_gc.consensus_stats"
        "theiacov_gc.auspice_json"
        "theiacov_gc.aligned_bai"
        "theiacov_gc.json_summary"

    The merged summary is under: "theiacov_gc.summaries_tsv" and "theiacov_gc.summaries_json"
    Output files should start with the sample name (sample01.file or sample01_file)
    """
    mkdir(f"{args.outdir}")
    theiacov_results = None
    theiacov_json = None
    for key, outputs in metadata["outputs"].items():
        task_name = key.replace("theiacov_gc.", "")
        if args.debug:
            print(f"Working on {task_name} outputs", file=sys.stderr)
        if task_name == "summaries_tsv":
            if args.debug:
                print(f"Copying {outputs} to {args.outdir}/complete-theiacov-results.tsv", file=sys.stderr)
            copy2(outputs, f"{args.outdir}/theiacov-results.tsv")
            theiacov_results = read_theiacov_results(outputs)
        elif task_name == "summaries_json":
            if args.debug:
                print(f"Copying {outputs} to {args.outdir}/complete-theiacov-results.json", file=sys.stderr)
            copy2(outputs, f"{args.outdir}/theiacov-results.json")
            theiacov_json = read_theiacov_results(outputs, is_json=True)
        else:
            for output in outputs:
                samplename = os.path.basename(output)
                
                for ext in OUTPUTS[task_name]['ext']:
                    samplename = samplename.replace(ext, "")

                if task_name == 'reads_dehosted':
                    samplename = samplename.replace("_R1", "")
                    samplename = samplename.replace("_R2", "")

                if samplename not in samples:
                    print(f"Unable to associate {output} with a sample", file=sys.stderr)
                    continue
                
                copy_path = f"{args.outdir}/{OUTPUTS[task_name]['folder']}"
                mkdir(copy_path)

                if args.debug:
                    print(f"Copying {output} to {copy_path}", file=sys.stderr)

                if task_name == 'kraken_report_dehosted':
                    copy2(output, f'{copy_path}/{samplename}_dehosted_kraken2_report.txt')
                else:
                    copy2(output, copy_path)
