# Running Trinity

[Trinity](http://trinityrnaseq.github.io/) is a _de novo_ assembler of RNA-seq reads which can handle pretty large datasets written by [Brian Haas](https://personal.broadinstitute.org/bhaas/) and collaborators. It also has a setting called Genome Guided Trinity which works by first aligning reads to the genome and binning the genome into segements and assembly reads only in a region - this helps cut down on paralogous genes being co-assembled and also can be more accurate and a little faster since the problem is now subdivided into smaller portions for the assembler.

Nearly all the info you need for running Trinity is on the website supported by the developed so it will not be repeated here, but instead will emphasize some points of clarification and

# Getting data ready for assembly

## FASTQ from sequencing center

These data should work just fine out of the gate. Remembering there are some which will be paired end and some single end.

## FASTQ from NCBI / SRA

Fix read names ...

# Knowing if your data are strand-specific RNAseq

One critical aspect is knowing if the data are strand specific and the organization of the reads as RF or FR. Guessing this [can be done](https://www.biostars.org/p/98756/) but it requires a reference genome. However the [Trinity documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Examine-Strand-Specificity) discuss this so it may require you to run a whole assembly iteration with Trinity and then return and re-map reads against this transcript assembly.

Several tools will help guess this, though generally this is only going to be guessable if you also have a genome to align the reads to.
- [RSeQC](http://rseqc.sourceforge.net/) - RNA-seq Quality Control Package has a tool called [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)

# Running Trinity on HPCC

For large read set (eg 200M reads or more) you will want to assign a lot of memory. The `intel` queue max memory that can be requested is 500 Gb and the `highmem` queue can request up to 1Tb (1000Gb).

# Inferring proteins from Trinity assembly

The tool [Transdecoder](https://transdecoder.github.io/) (also written by [Brian Haas](https://personal.broadinstitute.org/bhaas/)) can be used to infer open reading frames (ORFs) from transcript assembly.

```bash
module load transdecoder
TransDecoder.LongOrfs -t Trinity.fasta
```
