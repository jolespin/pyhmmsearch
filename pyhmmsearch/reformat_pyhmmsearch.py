#!/usr/bin/env python
import sys, os, argparse, gzip
from collections import defaultdict
from tqdm import tqdm
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.4.25"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.tsv> -o <output.tsv> -f <format> ".format(__program__)
    epilog = "Copyright 2024 New Atlantis Labs (jolespin@newatlantis.io)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/input.tsv. The results of `pyhmmsearch.py` [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-f", "--format", type=str, choices={"table", "pickle"}, help="Output format")
    parser.add_argument("--no_header",action="store_true", help="Input does not have header")
    parser.add_argument("-b", "--best_hits_only",action="store_true", help="Best hits only")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Input
    if opts.input == "stdin":
        f_input = sys.stdin 
    else:
        if opts.input.endswith(".gz"):
            f_input = gzip.open(opts.input, "rt")
        else:
            f_input = open(opts.input, "r")

    if not opts.no_header:
        next(f_input)
    
    if opts.best_hits_only:
        output = defaultdict(dict)
        for line in tqdm(f_input, desc="Reading PyHMMSearch"):
            line = line.strip()
            if line:
                id_protein, id_hmm, threshold, score, bias, best_domain_score, best_domain_bias, evalue = line.split("\t")
                score = float(score)
                evalue = float(evalue)
                update = True
                if id_protein in output:
                    existing_score = output[id_protein]["score"]
                    if score <= existing_score:
                        update = False
                if update:
                    output[id_protein]["id"] = id_hmm
                    # output[id_protein]["name"] = definition
                    output[id_protein]["evalue"] = evalue
                    output[id_protein]["score"] = score
        df_output = pd.DataFrame(output).T
    else:
        output = defaultdict(lambda: defaultdict(list))
        for line in tqdm(f_input, desc="Reading PyHMMSearch"):
            line = line.strip()
            if line:
                id_protein, id_hmm, threshold, score, bias, best_domain_score, best_domain_bias, evalue = line.split("\t")
                output[id_protein]["ids"].append(id_hmm)
                # output[id_protein]["names"].append(definition)
                output[id_protein]["evalues"].append(float(evalue))
                output[id_protein]["scores"].append(float(score))
        df_output = pd.DataFrame(output).T
        df_output.insert(0, "number_of_hits", df_output["ids"].map(len))
    df_output.index.name = "id_protein"

    if opts.format == "pickle":
        df_output.to_pickle(opts.output)
    else:
        df_output.to_csv(opts.output, sep="\t")

    if opts.input != "stdin":
        f_input.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
