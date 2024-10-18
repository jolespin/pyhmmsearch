#### Daily Change Log: 
* [2024.10.18] - Uses entry points for executables instead of copying scripts to bin/
* [2024.7.29] - Added `IsADirectoryError` if `--hmm_database` is a directory in `pyhmmsearch.py`
* [2024.7.29] - Fixed file extension error where an extra "." was added during glob when directory of HMMs was provided for `--hmm_database` in `serialize_hmm_models.py`
#### Pending: 
* Need to add support for length cutoff [BUSCO Issue #74](https://gitlab.com/ezlab/busco/-/issues/740)