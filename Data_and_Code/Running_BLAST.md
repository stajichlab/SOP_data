# How to Run BLAST on the Commandline

## Goals

Sounds like you're interested in running Blast on the command line. There are some really easy steps to complete this.

Have you run NCBI Blast using their website? [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), here you take your gene of interest and can blast it through a larger database, in this case all of NCBI.

## Preparing data

First you need to make sure you have all the material you need in the folder you want to work in. 

Create the directory you want to keep everything contained in and follow the next steps.

1. Download your sequence of interest. Is it DNA or protein? Is it from a specific gene? This will be a very important consideration for later.
* You first need to find this gene/protein through [NCBI](https://www.ncbi.nlm.nih.gov) or [uniprot](https://www.uniprot.org/). You may also have defined genes based on gene names which can be found often in organism-focused databases eg [SGD](https://yeastgenome.org), [FungiDB](https://fungidb.org/fungidb).
* You can download your gene/protein based on an accession direct from the cluster using [efetch](https://www.ncbi.nlm.nih.gov/books/NBK179288/) also see this [tutorial on eutils](https://astrobiomike.github.io/unix/ncbi_eutils).

``` 
 module load ncbi_edirect/16.8.20220329 
 efetch -db nuccore -id HM011534.1 -format fasta > PKS1.fna 
 efetch -db protein -id ADF36144.1 -format fasta > PKS1.faa 
```

* Now you need to download the single genome/proteome or database you want to search. Databases NCBI and see if there are protein or nucleotide files of your genome available. 
> An example of where we can find these files is [NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome/?term=Aspergillus+nidulans). Where we are looking up a genome for *Aspergillus nidulans*. 
> As you can see there are multiple options: genome, transcript, protein. You can right click and press copy address/link for the FastA formatted file of interest.
> Then in your command line type in `curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/205/GCF_000149205.2_ASM14920v2/GCF_000149205.2_ASM14920v2_protein.faa.gz` or curl -O whatever the hyperlink is on your end.
* You should have downloaded the above protein file.
> To obtain the nucleotide file, download the genome using the same steps above: `curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/205/GCF_000149205.2_ASM14920v2/GCF_000149205.2_ASM14920v2_genomic.fna.gz`
> These files are compressed with `gzip` - to use these files in the next steps you will need to uncompress these for for example: `gunzip GCF_000149205.2_ASM14920v2_genomic.fna.gz` and `gunzip GCF_000149205.2_ASM14920v2_protein.faa.gz`


## Now you have both a protein and genome types to work with. 
Depending on the data types (protein or dna) you want to search NCBI BLAST commands there are: blastn (DNA to DNA), blastP (protein to protein) blastx (DNA translated to protein), tblastn (protein to translated DNA), tblastx (translated DNA to translated DNA).  Familiarize yourself with these concepts before progressing.
You can complete same processes on the command line from the website [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) on the commandline using the UNIX module  on the HPCC `module load ncbi-blast`. If you are running on your laptop or your own [bioconda installed blast](https://bioconda.github.io/recipes/blast/README.html). 

Once loaded in you can check all the other programs by typing in -h next to each command, example.
`blastn -h` or `tblastn -h`

It will display a comprehensive manual for each command. If you can't find exactly how to work this, check online for resources.

## In order to properly run your BLAST search make sure you prepare your database.
In this case, your database will be the genome/protein file you have pulled from NCBI.

To create your database run, `makeblastdb -in YOUR_ORGANISM_GENOME.fna.fasta -dbtype nucl -out YOUR_ORGANISM` if you have a protein file change the -dbtype tag to `prot` instead.

Once this is done running (shouldn't take longer than 30 seconds) you are ready to blast through your genome/proteome.

Next line of code to run your nucleotide or protein BLAST will look like: (These are just some examples, read through the manual to fine-tune this code)

```
module load ncbi-blast/2.13.0+

blastn -db YOUR_ORGANISM -query YOUR_GENE.fna -evalue 1e-5 -outfmt 6 -out YOUR_ORGANISM.blastx.tab

blastp -db YOUR_ORGANISM -query YOUR_PROTEIN.faa -max_target_seqs 10 -evalue 1e-5 -num_threads 15 -outfmt '7 std qlen slen' -out YOUR_PROTEIN.tsv

```

Once you have your results, good luck!

-T
