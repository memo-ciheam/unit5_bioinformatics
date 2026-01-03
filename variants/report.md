# Session 5 – Mapping, variant calling and variant filtering

<details>
<summary><strong>🔴 Exercise 1 – FASTQ data inspection and quality assessment</strong></summary>

---

<details>
<summary><strong>Exercise 1.1 – Why are there two FASTQ files?</strong></summary>

### Question
Why are there two FASTQ files with the same name but different numeric suffixes?

### Answer
The two FASTQ files indicate **paired-end sequencing**.

- `SAMEA2569438.chr10_1.fastq.gz` contains **Read 1 (forward reads)**
- `SAMEA2569438.chr10_2.fastq.gz` contains **Read 2 (reverse reads)**

Each DNA fragment was sequenced from both ends, generating two corresponding reads.

</details>

---

<details>
<summary><strong>Exercise 1.2 – FASTQ format and differences between Read 1 and Read 2</strong></summary>

### Question
What are the differences between the two FASTQ files?  
Do they follow the typical FASTQ format?

### Answer
Both files follow the **standard FASTQ format** and show very similar characteristics:

- Same number of reads (**172,253**)
- Same GC content (**40%**)
- Similar read length distribution (**30–83 bp**)
- Same encoding (**Sanger / Illumina 1.9**)

The only difference is the **read orientation** (forward vs reverse).

**Read 1 – Basic statistics:**  
![FastQC basic statistics – Read 1](../images/ex5.1.png)

**Read 2 – Basic statistics:**  
![FastQC basic statistics – Read 2](../images/ex5.5.png)

</details>

---

<details>
<summary><strong>Exercise 1.3 – Sequencing error rate along reads</strong></summary>

### Question
Are read starts and read ends similar in terms of error rate?

### Answer
No. Sequencing quality is **higher at the beginning of reads** and **decreases toward the read ends**.

This pattern is observed in both Read 1 and Read 2 and is typical for Illumina sequencing data.

**Per-base quality – Read 1:**  
![Per-base quality – Read 1](../images/ex5.2.png)

**Per-base quality – Read 2:**  
![Per-base quality – Read 2](../images/ex5.6.png)

</details>

---

### Additional quality metrics

**Per-sequence quality scores – Read 1:**  
![Per-sequence quality – Read 1](../images/ex5.3.png)

**Read length distribution – Read 1:**  
![Read length distribution – Read 1](../images/ex5.4.png)

**Per-sequence quality scores – Read 2:**  
![Per-sequence quality – Read 2](../images/ex5.7.png)

**Read length distribution – Read 2:**  
![Read length distribution – Read 2](../images/ex5.8.png)


---

### Interpretation

The FASTQ files correspond to high-quality paired-end Illumina sequencing data.  
The observed quality decay toward read ends is expected and does not compromise downstream mapping and variant calling.

### References

FastQC Documentation – Quality control for high throughput sequence data
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Illumina Sequencing Technology – Sequencing by synthesis overview
https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html

OpenAI ChatGPT – used for language refinement and conceptual clarification of FASTQ format and FastQC quality assessment
Prompt-based assistance during Session 5, Exercise 1

</details>

