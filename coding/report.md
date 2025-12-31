
# Session 1_Annotation of coding sequences

<details>
<summary><strong>🔴 Exercise 1</strong></summary> 
 
## Question  

How many sequences have been formatted 
and how does this affect the E-value of BLAST searches?

## Objective

The aim of this exercise was to prepare a protein sequence collection  
for BLAST searches.

In addition, the objective was to understand how the size of the database  
influences the statistical significance of BLAST results,  
with particular emphasis on the E-value.

## Commands used  

```bash
(base) mehmetsonmez@MEMO:~$ docker start -ai bioinfo_iamz

vep@5dc71ff4216c:/home/vep$ cd /home/vep/test_data

vep@5dc71ff4216c:/home/vep/test_data$ ls

vep@5dc71ff4216c:/home/vep/test_data$ cp uniprot_Atha.fasta.gz /data/

vep@5dc71ff4216c:/home/vep/test_data$ cd /data

vep@5dc71ff4216c:/data$ ls

vep@5dc71ff4216c:/data$ gunzip uniprot_Atha.fasta.gz

vep@5dc71ff4216c:/data$ ls

vep@5dc71ff4216c:/data$ /home/vep/get_homologues/bin/ncbi-blast-2.16.0+/bin/makeblastdb \
-dbtype prot \
-in uniprot_Atha.fasta

```

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

No difficulties

## References

BLAST Help Manual – E-value definition
https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html

OpenAI ChatGPT – used for language refinement and conceptual 
clarification of BLAST E-value interpretation
Prompt-based assistance during Session 1, Exercise 1

</details>

<details>
<summary><strong>🔴 Exercise 2</strong></summary>

## Question  

Can you redirect the output of `blastp` and `blastx`  
to separate files called `test.faa.blast` and `test.fna.blast`?

## Objective  

The objective of this exercise was to query the UniProt protein collection  
using a sample coding sequence (ARF6) and its corresponding transcript.

In addition, the aim was to learn how to redirect BLAST search outputs  
from the standard output to separate result files for further analysis.

## Commands used  

```bash
(base) mehmetsonmez@MEMO:~$ docker start -ai bioinfo_iamz

vep@5dc71ff4216c:/home/vep$ cd /home/vep/test_data

vep@5dc71ff4216c:/home/vep/test_data$ cp test.faa /data/

vep@5dc71ff4216c:/home/vep/test_data$ cp test.fna /data/

vep@5dc71ff4216c:/home/vep/test_data$ cd /data

vep@5dc71ff4216c:/data$ /home/vep/get_homologues/bin/ncbi-blast-2.16.0+/bin/blastp \
-db uniprot_Atha.fasta \
-query test.faa \
-outfmt 6 \
> test.faa.blast

vep@5dc71ff4216c:/data$ /home/vep/get_homologues/bin/ncbi-blast-2.16.0+/bin/blastx \
-db uniprot_Atha.fasta \
-query test.fna \
-outfmt 6 \
> test.fna.blast

vep@5dc71ff4216c:/data$ ls
```
## Results  

Two separate BLAST output files were successfully generated:

- `test.faa.blast`, containing the results of the protein–protein BLAST search (`blastp`)
- `test.fna.blast`, containing the results of the transcript–protein BLAST search (`blastx`)

Both files are in tab-delimited format (`outfmt 6`), which is suitable  
for downstream filtering and analysis.

## Interpretation and discussion  

The redirection operator (`>`) allows BLAST results to be written directly  
to files instead of being printed to the terminal.

This approach is essential when working with large BLAST outputs,  
as it enables reproducible analyses, result storage, and further processing  
without losing information displayed on the screen.

Using separate output files for `blastp` and `blastx` ensures that  
protein-based and transcript-based searches are clearly distinguished.

## Difficulties encountered  

The BLAST executables were not available in the default system PATH  
inside the Docker container. Therefore, the full path to the BLAST binaries  
had to be used in order to execute the commands successfully.

## References  

BLAST Help Manual  
https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html

OpenAI ChatGPT – used for language refinement and clarification  
of BLAST output redirection during Session 1, Exercise 2

</details>
