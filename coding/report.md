
# Session 1_Annotation of coding sequences

<details>
<summary><strong>Exercise 1</strong></summary> 
 
## Question  
**How many sequences have been formatted 
and how does this affect the E-value of BLAST searches?**

## Objective
The aim of this exercise was to prepare a protein sequence collection  
for BLAST searches.

In addition, the objective was to understand how the size of the database  
influences the statistical significance of BLAST results,  
with particular emphasis on the E-value.

## Commands used  

First, the compressed UniProt protein dataset was decompressed:

```bash
gunzip uniprot_Atha.fasta.gz
```
Then, the protein FASTA file was formatted as a BLAST database 
using the makeblastdb command:
```bash
/home/vep/get_homologues/bin/ncbi-blast-2.16.0+/bin/makeblastdb \
  -dbtype prot \
  -in uniprot_Atha.fasta
```
This command creates all the index files required by BLAST to
efficiently search protein sequences.

## Results

The database formatting process produced the following output:

Adding sequences from FASTA; added **15719 sequences in 0.166255 seconds**.

This indicates that a total of 15,719 protein sequences from 
Arabidopsis thaliana were successfully included in the BLAST database.


## Interpretation and discussion

The E-value in BLAST describes the number of alignments with a given 
score that are expected to occur by chance when searching a database 
of a specific size. Therefore, the E-value is directly dependent on 
the size of the database.

In this exercise, formatting a database containing 15,719 protein 
sequences means that BLAST searches are performed against a relatively
large search space. As a consequence, for a given alignment score, 
the E-value will be higher than it would be in a smaller database. 
This is because the probability of finding random matches increases 
as the number of sequences in the database increases.

Thus, when interpreting BLAST results, it is essential to consider 
database size: a low E-value in a large database provides stronger
evidence of biological relevance than the same score obtained from a 
small database.

## Difficulties encountered

The main difficulty encountered during this exercise was related to the
execution environment. The BLAST software was installed inside a Docker
container, but the makeblastdb executable was not available in the 
default system PATH. As a result, the correct absolute path to the 
BLAST binaries had to be identified before the command could be executed
successfully.

## References

BLAST Help Manual – E-value definition
https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html

(visited 15/12/2025)

OpenAI ChatGPT – used for language refinement and conceptual 
clarification of BLAST E-value interpretation
Prompt-based assistance during Session 1, Exercise 1

</details>

