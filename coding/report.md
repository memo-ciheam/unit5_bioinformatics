## Exe1 — Formatting a sequence collection for BLAST

The following command was used to format the protein database:

```bash
makeblastdb -dbtype prot -in uniprot_Atha.fasta

The output indicates that 15,719 protein sequences were added to the BLAST database.

According to the BLAST documentation, the E-value describes the number of hits expected by chance when searching a database of a given size. Therefore, as the database size increases, the E-value for a given alignment score also increases, because the probability of finding random matches becomes higher.

Reference
https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html



