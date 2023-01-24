makeblastdb -in transgenicRef.fasta -dbtype nucl -parse_seqids 
blastn -query queryprimers.fasta -subject transgenicRef.fasta -evalue 1e-23 -outfmt 6 -out primerblast.txt -max_target_seqs 2