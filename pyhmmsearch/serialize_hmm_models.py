#!/usr/bin/env python
import sys, os, glob, gzip, warnings, argparse, pickle
from collections import defaultdict
from tqdm import tqdm
from pyhmmer.plan7 import HMMFile

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.7.29"

def stdin_is_empty():
    return sys.stdin.isatty()

def main(args=None):
    # Options
    # =======
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -d <path/to/hmm_database> -b <hmm_database.pkl.gz> -f name ".format(__program__)
    epilog = "Copyright 2024 New Atlantis Labs (jolespin@newatlantis.io)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_database = parser.add_argument_group('Database arguments')
    parser_database.add_argument("-l", "--hmm_list", default="stdin", type=str, help="path/to/hmm_database.  If it is a single file then it is interpreted as a multi-HMM database. If it's a directory, then it will be scanned for .hmm[.gz] files (Not case sensitive). Cannot be used with --hmm_database.")
    parser_database.add_argument("-d", "--hmm_database", type=str, help="path/to/hmm_database.  If it is a single file then it is interpreted as a multi-HMM database. If it's a directory, then it will be scanned for .hmm[.gz] files (Not case sensitive). Cannot be used with --hmm_list.")
    parser_database.add_argument("-b", "--serialized_database", required=True, type=str, help="path/to/database.pkl[.gz] is dictionary of HMM models")
    parser_database.add_argument("-f", "--hmm_marker_field", default="accession", type=str, help="HMM reference type (accession, name) [Default: accession]")

    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Input
    # =====
    hmm_filepaths = list()
    if opts.hmm_database:
        assert opts.hmm_list == "stdin", "Cannot provide --hmm_database and --hmm_list"
        assert stdin_is_empty(), "Cannot provide stdin if using --hmm_database"
        if os.path.isdir(opts.hmm_database):
            for extension in ["hmm", "hmm.gz", "HMM", "HMM.gz"]:
                for fp in glob.glob(os.path.join(opts.hmm_database, f"*.{extension}")):
                    hmm_filepaths.append(fp)
        else:
            hmm_filepaths.append(opts.hmm_database)
    else:
        if opts.hmm_list == "stdin":
            f_in = sys.stdin
        else:
            if opts.hmm_list.endswith(".gz"):
                f_in = gzip.open(opts.hmm_list, "rt")
            else:
                f_in = open(opts.hmm_list, "r")
            
        for line in f_in:
            fp = line.strip()
            if fp:
                hmm_filepaths.append(fp)
        if f_in != sys.stdin:
            f_in.close()
            

    # Database
    # ========
    # Load HMMs
    name_to_hmm = dict()
    for fp in hmm_filepaths:
        with HMMFile(fp) as f_hmm:
            hmms = list(f_hmm)
            for hmm in tqdm(hmms, desc=f"Loading HMMs from {fp}", total=len(hmms)):
                try:
                    name = getattr(hmm, opts.hmm_marker_field)
                except AttributeError:
                    if opts.hmm_marker_field == "name":
                        warnings.warn("HMM did not have `--hmm_marker_field name`.  Try using `--hmm_marker_field accession")
                    if opts.hmm_marker_field == "accession":
                        warnings.warn("HMM did not have `--hmm_marker_field accession`.  Try using `--hmm_marker_field name`")
                    sys.exit(1)
                name = name.decode()
                name_to_hmm[name] = hmm

    # Output
    # ======
    # Write serialized database
    if opts.serialized_database.endswith((".gz", ".pgz")):
        f_out = gzip.open(opts.serialized_database, "wb")
    else:
        f_out = open(opts.serialized_database, "wb")
    pickle.dump(name_to_hmm, f_out)
    f_out.close()


if __name__ == "__main__":
    main(sys.argv[1:])
