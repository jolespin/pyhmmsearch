# PyHMMSearch
Fast implementation of [`HMMSearch`](https://github.com/EddyRivasLab/hmmer) optimized for high-memory systems using [`PyHmmer`](https://github.com/althonos/pyhmmer).  `PyHMMSearch` can handle fasta in uncompressed or gzip format and databases in either HMM or Python pickle serialized format.  No intermediate files are created. 

### Installation:

```
pip install pyhmmsearch
``` 

### Dependencies: 

* pyhmmer >=0.10.12
* pandas
* tqdm


### Benchmarking:

| Database | Tool            | Single Threaded | 12 Threads |
|----------|-----------------|-----------------|------------|
| Pfam     | PyHMMSearch     | 2:24            | 0:20       |
| Pfam     | HMMER HMMSearch | 2:53            | 2:27       |

\* Time in minutes for 4977 proteins in `test/test.faa.gz`.  

Official benchmarking for `hmmsearch` algorithm implemented in `PyHMMER`  against `HMMER` from [Larralde et al. 2023](https://academic.oup.com/bioinformatics/article/39/5/btad214/7131068): 


<img src="images/pyhmmer-benchmarking.jpeg" alt="drawing" width="500"/>



### Usage:
Recommended usage for `PyHMMSearch` is on systems with 1) high RAM;  2) large numbers of threads; and/or 3) reading/writing to disk is charged (e.g., AWS EFS).  Also useful when querying a large number of proteins. 

* #### Using the official Pfam database files:

    ```bash
    # Download database
    DATABASE_DIRECTORY=/path/to/database_directory/
    mkdir -p ${DATABASE_DIRECTORY}/Annotate/Pfam
    wget -v -P ${DATABASE_DIRECTORY}/Annotate/Pfam https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 

    # Run PyHMMSearch
    pyhmmsearch -i test/test.faa.gz  -o output.tsv -b ${DATABASE_DIRECTORY}/Annotate/Pfam/Pfam-A.hmm.gz -p=-1
    ```

* #### Build a serialized database:

    ```bash
    # Provide a database
    serialize_hmm_models -d path/to/Pfam-A.hmm.gz  -b path/to/database.pkl.gz

    # or a directory of HMMs
    serialize_hmm_models -d path/to/hmm_directory/  -b path/to/database.pkl.gz

    # or from a list of filepaths to HMM models
    serialize_hmm_models -l path/to/hmms.list  -b path/to/database.pkl.gz

    # or form a list through stdin
    ls path/to/directory/*.hmm |  serialize_hmm_models -b path/to/database.pkl.gz
    ```

* #### Using the serialized database files:

    Database can be uncompressed pickle or gzipped pickle. 

    ```bash
    pyhmmsearch -i test/test.faa.gz  -o output.tsv -b ~/Databases/Pfam/database.pkl.gz -p=-1
    ```

* #### Using a custom threshold table (e.g., BUSCO Markers):

    ```bash
    pyhmmsearch -i test/test.faa.gz  -o output.tsv -d test/bacteria_odb10/bacteria_odb10.hmm.gz -s test/bacteria_odb10/scores_cutoff -f name -p=-1
    ```

* #### Grouping hits by query protein:

    ```bash
    reformat_pyhmmsearch -i pyhmmsearch_output.tsv -o pyhmmsearch_output.reformatted.tsv
    ```

#### Options:

```
$ pyhmmsearch -h
usage: pyhmmsearch -i <proteins.fasta> -o <output.tsv> -d 

    Running: pyhmmsearch v2025.10.23 via Python v3.14.0 | /Users/jolespin/miniforge3/envs/annotate/bin/python3.14

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

I/O arguments:
  -i, --proteins PROTEINS
                        path/to/proteins.fasta. stdin does not stream and loads everything into memory. [Default: stdin]
  -o, --output OUTPUT   path/to/output.tsv [Default: stdout]
  --no_header           No header
  --tblout TBLOUT       path/to/output.tblout
  --domtblout DOMTBLOUT
                        path/to/output.domtblout

Utility arguments:
  -p, --n_jobs N_JOBS   Number of threads to use [Default: 1]

HMMSearch arguments:
  -s, --scores_cutoff SCORES_CUTOFF
                        path/to/scores_cutoff.tsv[.gz] [id_hmm]<tab>[score_threshold], No header.
  -f, --hmm_marker_field {accession,name}
                        HMM reference type (accession, name) [Default: accession]
  -t, --score_type {domain,full}
                        {full, domain} [Default: full]
  -m, --threshold_method {gathering,e,noise,trusted}
                        Cutoff threshold method [Default:  e]
  -e, --evalue EVALUE   E-value threshold [Default: 10.0]

Database arguments:
  -d, --hmm_database HMM_DATABASE
                        path/to/database.hmm cannot be used with -b/-serialized_database.  Expects a (concatenated) HMM file and not a directory. You can build a database from a directory using `serialize_hmm_models.py`
  -b, --serialized_database SERIALIZED_DATABASE
                        path/to/database.pkl cannot be used with -d/--database_directory.  Database should be pickled dictionary {name:hmm}

https://github.com/jolespin/pyhmmsearch
```


#### Outputs:

* From **pyhmmsearch**:

    | id_protein | id_hmm | threshold | score | bias | best_domain-score | best_domain-bias | e-value | description |
    |------------|--------|-----------|-------|------|-------------------|------------------|---------|-------------|
    | SRR13615825__k127_409025_2 | PF09847.13 |  | 13.46 | 0.609 | 13.226 | 0.609 | 0.008249 | Membrane protein of 12 TMs |
    | SRR13615825__k127_453760_1 | PF00389.34 |  | 93.686 | 6.702 | 89.856 | 6.702 | 1.9839999999999998e-27 | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain |
    | SRR13615825__k127_295655_1 | PF00389.34 |  | 83.195 | 0.005 | 83.167 | 0.005 | 3.456e-24 | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain |
    | SRR13615825__k127_218710_3 | PF00389.34 |  | 42.235 | 0.004 | 42.073 | 0.004 | 1.559e-11 | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain |
    | SRR13615825__k127_272080_1 | PF00389.34 |  | 24.673 | 0.0 | 22.067 | 0.0 | 4.154e-06 | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain |

* From **reformat_pyhmmsearch**:

    | id_protein                 | number_of_hits | ids                                                      | evalues                                     | scores                            | descriptions                                                                                                                                                                                                                         |
    |----------------------------|----------------|----------------------------------------------------------|---------------------------------------------|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | SRR13615825__k127_409025_2 | 2              | ['PF09847.13', 'PF02659.19']                             | [0.008249, 1.763e-13]                       | [13.46, 48.757]                   | ['Membrane protein of 12 TMs', 'Putative manganese efflux pump']                                                                                                                                                                     |
    | SRR13615825__k127_453760_1 | 3              | ['PF00389.34', 'PF02826.23', 'PF03446.19']               | [1.984e-27, 2.113e-39, 2.41e-08]            | [93.686, 132.902, 32.336]         | ['D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain', 'D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain', 'NAD binding domain of 6-phosphogluconate dehydrogenase']                                   |
    | SRR13615825__k127_295655_1 | 2              | ['PF00389.34', 'PF02826.23']                             | [3.456e-24, 7.794e-21]                      | [83.195, 72.421]                  | ['D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain', 'D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain']                                                                                             |
    | SRR13615825__k127_218710_3 | 1              | ['PF00389.34']                                           | [1.559e-11]                                 | [42.235]                          | ['D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain']                                                                                                                                                                  |
    | SRR13615825__k127_272080_1 | 4              | ['PF00389.34', 'PF02826.23', 'PF03720.19', 'PF13478.10'] | [4.154e-06, 2.035e-41, 0.005561, 0.0005575] | [24.673, 139.471, 15.248, 18.593] | ['D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain', 'D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain', 'UDP-glucose/GDP-mannose dehydrogenase family, UDP binding domain', 'XdhC Rossmann domain'] |

* From **reformat_pyhmmsearch** with **-b/--best_hits_only**:

    | id_protein                 | id         | evalue                 | score   | description                                                       |
    |----------------------------|------------|------------------------|---------|-------------------------------------------------------------------|
    | SRR13615825__k127_409025_2 | PF02659.19 | 1.763e-13              | 48.757  | Putative manganese efflux pump                                    |
    | SRR13615825__k127_453760_1 | PF02826.23 | 2.113e-39              | 132.902 | D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain |
    | SRR13615825__k127_295655_1 | PF00389.34 | 3.456e-24              | 83.195  | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain   |
    | SRR13615825__k127_218710_3 | PF00389.34 | 1.559e-11              | 42.235  | D-isomer specific 2-hydroxyacid dehydrogenase, catalytic domain   |
    | SRR13615825__k127_272080_1 | PF02826.23 | 2.0349999999999998e-41 | 139.471 | D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain |

#### If you use this tool, please cite the following sources: 

* Eddy SR. Accelerated Profile HMM Searches. PLoS Comput Biol. 2011 Oct;7(10):e1002195. doi: 10.1371/journal.pcbi.1002195. Epub 2011 Oct 20. PMID: 22039361; PMCID: PMC3197634.

* Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics. 2023 May 4;39(5):btad214. doi: 10.1093/bioinformatics/btad214. PMID: 37074928; PMCID: PMC10159651.


#### License:

The code for `PyHMMSearch` is licensed under an MIT License

Please contact jol.espinoz@gmail.com regarding any licensing concerns.
