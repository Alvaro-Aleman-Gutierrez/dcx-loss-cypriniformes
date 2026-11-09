#!/bin/bash
# Sintenia del bloque DCX en lochas: Nemacheilidae con DCX (Barbatula, Triplophysa) vs Cobitidae sin DCX (Misgurnus, Paramisgurnus).
P=${PROJECT}; O=${DATA}/syn_loach; B=${DATA}/p1b/genomes
cd $O
for sp in Barbatula_barbatula Triplophysa_rosa Misgurnus_anguillicaudatus Paramisgurnus_dabryanus; do
  db=$P/genomes/db_$sp; [ -e $db.nsq ] || [ -e $db.ndb ] || db=$B/db_$sp
  [ -s $sp.blk.tsv ] || tblastn -query $P/syn/block_q.faa -db $db -evalue 1e-10 -outfmt "6 qseqid sseqid pident length sstart send bitscore" \
     -num_threads 2 -max_target_seqs 5 -max_hsps 3 > $sp.blk.tsv
  [ -s $sp.dcx.tsv ] || tblastn -query $P/tbn/q_dcx.faa -db $db -evalue 1e-5 -outfmt "6 qseqid sseqid pident length qstart qend sstart send bitscore sseq" \
     -num_threads 2 -max_target_seqs 50 > $sp.dcx.tsv
done
echo FIN > done.txt
