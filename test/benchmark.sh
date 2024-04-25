echo "PyHMMSearch 12 threads"
time python ../pyhmmsearch/pyhmmsearch.py -i test.faa -p 12 -o pyhmmsearch_cpu12.tsv -b /Users/jolespin/Databases/Pfam/database.pkl.gz -m gathering
echo "PyHMMSearch 1 thread"
time python ../pyhmmsearch/pyhmmsearch.py -i test.faa -p 1 -o pyhmmsearch.tsv -b /Users/jolespin/Databases/Pfam/database.pkl.gz -m gathering

echo "HMMSearch 12 threads"
time hmmsearch -o /dev/null --tblout hmmsearch_cpu12.tsv --cut_ga --cpu 12 /Users/jolespin/Databases/Pfam/Pfam-A.hmm.gz test.faa
echo "HMMSearch 1 thread"
time hmmsearch -o /dev/null --tblout hmmsearch.tsv --cut_ga --cpu 1 /Users/jolespin/Databases/Pfam/Pfam-A.hmm.gz test.faa
