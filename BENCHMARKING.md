## PyHMMSearch

### Pfam Database

#### 12 threads

```
$ time python ../pyhmmsearch/pyhmmsearch.py -i test.faa -p 12 -o pyhmmsearch_cpu12.tsv -b /Users/jolespin/Databases/Pfam/database.pkl.gz -m gathering

Loading serialized HMM database
Performing HMMSearch: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20795/20795 [00:16<00:00, 1297.51it/s]

real	0m20.667s
user	3m0.967s
sys	0m2.842s
```

#### Single thread
```
$ time python ../pyhmmsearch/pyhmmsearch.py -i test.faa -p 1 -o pyhmmsearch.tsv -b /Users/jolespin/Databases/Pfam/database.pkl.gz -m gathering

Loading serialized HMM database
Performing HMMSearch: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20795/20795 [02:20<00:00, 148.48it/s]

real	2m24.450s
user	2m20.839s
sys	0m3.328s
```

_________________________________________________________
## HMMER's HMMSearch

### Pfam Database

#### 12 threads

```
$ time hmmsearch -o /dev/null --tblout hmmsearch_cpu12.tsv --cut_ga --cpu 12 /Users/jolespin/Databases/Pfam/Pfam-A.hmm.gz test.faa

real	2m26.594s
user	4m7.765s
sys	0m51.786s
```

#### Single thread

```
$ time hmmsearch -o /dev/null --tblout hmmsearch.tsv --cut_ga --cpu 1 /Users/jolespin/Databases/Pfam/Pfam-A.hmm.gz test.faa

real	2m53.103s
user	3m54.172s
sys	0m8.797s
```

