#### Daily Change Log: 
* [2025.10.23] - Added `description` field to output of `pyhmmsearch` and `descriptions` to `reformat_pyhmmsearch`
* [2025.9.5] - Pinned `pyhmmer` version >=0.10.2,<0.11
* [2025.9.4] - Added `biopython` to `requirements.txt`
* [2025.1.23] - Added `--tblout` and `--domtblout` output arguments to main executable
* [2024.10.18] - Uses entry points for executables instead of copying scripts to bin/
* [2024.7.29] - Added `IsADirectoryError` if `--hmm_database` is a directory in `pyhmmsearch.py`
* [2024.7.29] - Fixed file extension error where an extra "." was added during glob when directory of HMMs was provided for `--hmm_database` in `serialize_hmm_models.py`

#### Pending: 
* Need to add support for length cutoff [BUSCO Issue #74](https://gitlab.com/ezlab/busco/-/issues/740)
* Add alignment coverage to output
