#!/usr/bin/env python
import sys, os, glob, gzip, warnings, argparse, pickle
from collections import defaultdict
from multiprocessing import cpu_count
from tqdm import tqdm
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, TextSequence, Alphabet
from pyhmmer import hmmsearch

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.4.25"

# Filter 
def filter_hmmsearch_threshold(
    hit, 
    threshold,
    score_type, 
    ):
    if score_type:
        if score_type == "domain":
            score = hit.best_domain.score
        else:
            score = hit.score
        if score >= threshold:
            evalue = hit.evalue
            return (threshold, score, evalue)
    # else:
    #     score = hit.score
    #     evalue = hit.evalue
    #     return ("", score, evalue)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <proteins.fasta> -o <output.tsv> -d ".format(__program__)
    epilog = "Copyright 2024 Josh L. Espinoza (jolespin@newatlantis.io)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i","--proteins", type=str, default="stdin", help = "path/to/proteins.fasta. stdin does not stream and loads everything into memory. [Default: stdin]")
    parser_io.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser_io.add_argument("--no_header", action="store_true", help = "No header")

    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-p","--n_jobs", type=int, default=1,  help = "Number of threads to use [Default: 1]")

    parser_hmmsearch = parser.add_argument_group('HMMSearch arguments')
    parser_hmmsearch.add_argument("-s", "--scores_cutoff", type=str, help="path/to/scores_cutoff.tsv [id_hmm]<tab>[score_threshold], No header.")
    parser_hmmsearch.add_argument("-f", "--hmm_marker_field", default="accession", type=str, choices={"accession", "name"}, help="HMM reference type (accession, name) [Default: accession]")
    parser_hmmsearch.add_argument("-t", "--score_type",  default="full", type=str, help="{full, domain} [Default: full]")
    parser_hmmsearch.add_argument("-m", "--threshold_method", type=str, default="e",choices={"gathering", "noise", "trusted", "e"},  help="Cutoff threshold method [Default:  e]")
    parser_hmmsearch.add_argument("-e","--evalue", type=float, default=10.0,  help = "E-value threshold [Default: 10.0]")

    parser_database = parser.add_argument_group('Database arguments')
    parser_database.add_argument("-d", "--hmm_database", type=str, help="path/to/database.hmm cannot be used with -b/-serialized_database")
    parser_database.add_argument("-b", "--serialized_database", type=str, help="path/to/database.pkl cannot be used with -d/--database_directory.  Database should be pickled dictionary {name:hmm}")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Threads
    if opts.n_jobs < 0:
        opts.n_jobs = cpu_count()

    # Database
    if opts.serialized_database:
        print("Loading serialized HMM database", file=sys.stderr)
        # Load serialized database
        if opts.serialized_database.endswith((".gz", ".pgz")):
            f = gzip.open(opts.serialized_database, "rb")
        else:
            f = open(opts.serialized_database, "rb")
        name_to_hmm = pickle.load(f)
        f.close()

    else:
                
        # Load HMMs
        name_to_hmm = dict()

        with HMMFile(opts.hmm_database) as f:
            for hmm in list(f):
                name = getattr(hmm, opts.hmm_marker_field)
                assert name is not None, "-f/--hmm_marker_field {}` returned NoneType for HMM.  Try `-f/--hmm_marker_field {}` instead".format(opts.hmm_marker_field, {"accession":"name", "name":"accession"}[opts.hmm_marker_field])
                name_to_hmm[name.decode()] = hmm

    if opts.scores_cutoff:
        name_to_threshold = dict()
        with open(opts.scores_cutoff, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    name, threshold = line.split("\t")
                    if name in name_to_hmm:
                        name_to_threshold[name] = float(threshold)
        A = set(name_to_hmm.keys())
        B = set(name_to_threshold.keys())
        assert A == B, "There are {} HMMs that do not have a score threshold".format(len(A - B))

    # Output
    if opts.output == "stdout":
        f_output = sys.stdout 
    else:
        f_output = open(opts.output, "w")

    if not opts.no_header:
        print("id_protein", "id_hmm", "threshold", "score", "bias", "best_domain-score", "best_domain-bias", "e-value", sep="\t", file=f_output)

    # Input
    if opts.proteins == "stdin":
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        proteins = list()
        for header, seq in tqdm(SimpleFastaParser(sys.stdin), f"Parsing sequences from {sys.stdin}"):
            id = header.split(" ")[0]
            digital_sequence = TextSequence(sequence=seq, name=id.encode()).digitize(Alphabet.amino())
            proteins.append(digital_sequence)

    else:
        with SequenceFile(opts.proteins, format="fasta", digital=True) as f:
            proteins = f.read_block()#sequences=opts.sequences_per_block)

    # Run HMMSearch  
    params = {
        "E":opts.evalue,
    }
    if opts.threshold_method != "e":
        params["bit_cutoffs"] = opts.threshold_method

    if opts.scores_cutoff:
        for id_hmm, hits in tqdm(zip(name_to_hmm.keys(), hmmsearch(name_to_hmm.values(), proteins, cpus=opts.n_jobs, **params)), desc="Performing HMMSearch", total=len(name_to_hmm)):
            for hit in hits:
                if hit.included:
                    result = filter_hmmsearch_threshold(
                        hit, 
                        threshold=name_to_threshold[id_hmm], 
                        score_type=opts.score_type,
                        )
                    if result:
                        id_query = hit.name.decode()

                        threshold, score, evalue = result
                        print(
                            id_query, 
                            id_hmm, 
                            threshold, 
                            "{:0.3f}".format(hit.score), 
                            "{:0.3f}".format(hit.bias), 
                            "{:0.3f}".format(hit.best_domain.score), 
                            "{:0.3f}".format(hit.best_domain.bias), 
                            "{:0.3e}".format(evalue), 
                        sep="\t", 
                        file=f_output,
                        )
    else:
        if opts.threshold_method == "e":
            for id_hmm, hits in tqdm(zip(name_to_hmm.keys(), hmmsearch(name_to_hmm.values(), proteins, cpus=opts.n_jobs, **params)), desc="Performing HMMSearch", total=len(name_to_hmm)):
                for hit in hits:
                    if hit.included:
                            id_query = hit.name.decode()
                            print(
                                id_query, 
                                id_hmm, 
                                "", 
                                "{:0.3f}".format(hit.score), 
                                "{:0.3f}".format(hit.bias), 
                                "{:0.3f}".format(hit.best_domain.score), 
                                "{:0.3f}".format(hit.best_domain.bias), 
                                "{:0.3e}".format(hit.evalue), 
                            sep="\t", 
                            file=f_output,
                            )
        else:
            for id_hmm, hits in tqdm(zip(name_to_hmm.keys(), hmmsearch(name_to_hmm.values(), proteins, cpus=opts.n_jobs, **params)), desc="Performing HMMSearch", total=len(name_to_hmm)):
                for hit in hits:
                    if hit.included:
                        id_query = hit.name.decode()
                        hmm = name_to_hmm[id_hmm]
                        threshold = getattr(hmm.cutoffs, opts.threshold_method)
                        print(
                            id_query, 
                            id_hmm, 
                            threshold, 
                            "{:0.3f}".format(hit.score), 
                            "{:0.3f}".format(hit.bias), 
                            "{:0.3f}".format(hit.best_domain.score), 
                            "{:0.3f}".format(hit.best_domain.bias), 
                            "{:0.3e}".format(hit.evalue), 
                        sep="\t", 
                        file=f_output,
                        )

    if f_output != sys.stdout:
        f_output.close()

if __name__ == "__main__":
    main(sys.argv[1:])
