
# Creating A BioProject

BioProjects are summary records for individual or multiple genome, amplicon, or transcriptome datasets. BioProjects can also be organized as umbrella projects so that a single BioProject refers to multiple other bioprojects. Data can be assocated with a Bioproject like Samples, SRA records, Genomes, and Publications. BioSamples are individual project and can be organized 1:1 for a single strain genome and a single Bioproject. Or as a multi-sample project which can be multiple samples from the same species (population genomics) and in fact multiple species. A single biosample can also be associated with multiple BioProjects where appropriate.

One important aspect is __Locus Prefix__. Make sure you have checked that box as we will use that as part of annotation steps. If you forget to check this you can request it by email from ncbi curators - see the email at the bottom of the bioproject page for your record.

Generally I don't see a problem making a BioProject public right away, but some prefer to wait until they are ready to release their data publically. Generally you can embargo genome and SRA records even if a BioProject is public.

![NCBI submit page 1](../img/ncbi_submit_bioproject_1.png)

# Creating BioSamples

A BioSample record is needed for each biological sample of material that was used. For example a sample to represent the strain where DNA was prepared from. It is not necessary to make a separate biosample if the same strain was grown several times and sequenced more than once (eg a Nanopore and Illumina run). It is typical for labs using the same strain to still deposit a biosample for each of these to reflect there may be some evolution of strains throughout.  However for RNA sequencing or any other control/treatment experiments each sample will need a separate biosample in order to associate the conditions of the sample.

# Depositing data in NCBI

To deposit raw data from sequencing projects in NCBI SRA or from annotation and assembly projects we stage the data in an upload from [HPCC](https://hpcc.ucr.edu). You can also do direct uploads via NCBI web interface which will use aspera if you choose, however the data must then be local to your desktop/laptop.  To avoid copying the data first to your computer can simply direct upload these data from the cluster/linux machines where your analyses were are performed.

To stage data with an upload to the NCBI upload site we use a command line tool. It can be done with regular FTP upload.
However for larger dataset using [Aspera](https://www.ncbi.nlm.nih.gov/books/NBK242625/) is better as it is a faster upload mechanism which utilizes a parallel networking technology called [FASP](https://en.wikipedia.org/wiki/Fast_and_Secure_Protocol). It can achieve upload speeds 10-50x faster than typical single FTP upload speeds so for large Illumina datasets this is critical.

To do an upload you need to use  out your unique folder associated with your NCBI account you need to go to [NCBI submission] (https://submit.ncbi.nlm.nih.gov/subs/genome/) site and click on either the Aspera command line or FTP menu.
![NCBI submit page 1](../img/ncbi_submit_main.png)

![NCBI submit page 2](../img/ncbi_submit_main_2.png)

```
module load aspera
# (all one line)
ascp -i $ASPERAKEY -QT -l500m -k1
-d <path/to/folder/containing files/>
subasp@upload.ncbi.nlm.nih.gov:uploads/YOUREMAIL_XXXX/PRJYYY
```

Alternatively to use FTP you can use ftp client.
```
# lftp is easy and can be put in the background
lftp subftp@ftp-private.ncbi.nlm.nih.gov
password: <put in your password>
> cd uploads/UNIQUE_PATH_FROM_NCBI_PAGE
> mkdir PRJYYY
> cd PRJYY
> mput *.fastq.gz
```

These data are pulled in by the NCBI in the upload process, but it will take 10-15 minutes once you have uploaded the data. If you have datasets that are associated with different upload projects (eg Genome 1, Genome 2, Amplicon data set 2) make a folder for each genome dataset because during the submission upload you will select a folder to associate with a submit. I usually make a folder based on each Bioproject - so the upload folder name can be in the upload name (since it will have the bioproject name) `subasp@upload.ncbi.nlm.nih.gov:uploads/YOUREMAIL_XXXX/BIOPROJECTXX`.

![NCBI submit Genome page 1](../img/ncbi_submit_genome_1.png)

![NCBI submit Genome page 2](../img/ncbi_submit_genome_2.png)

# SRA deposition

Amplicon data

## Raw DNA Sequence

Illumina data

PacBio and Nanopore

These are also similar to Illumina data deposition. Just FASTQ files can be desposited. However fast5 files are also deposited for Nanopore results.


## Raw RNA and Epigenomic data deposition

It is more typical to upload RNAseq data through the [Gene Expression Ombnibus](https://www.ncbi.nlm.nih.gov/geo/) (GEO). GEO will take uploaded raw data as well as computed gene expression data. It can also take ChIP-Seq Epigenomics

# Genome deposition

Genomes can be deposited as simply and assembly (FastA format) or an assembly with gene annotation (produced from a tool like funannotate). Annotated genomes in GenBank are usually deposited in the sequin format (`.sqn`). Typically a tool like funannotate will produce this already but it is generated from a simple annotation format called NCBI Table format (`.tbl`) and a genome assembly Fasta file. The `.tbl + .fsa` files are combined in a tool called `tbl2asn` which produces the `.sqn` file and does a bunch of checking for everything from frame shifts to ensuring gene/mRNA/exon features all overlap as expected in their ranges, and many more checks.

The first page for genome submission is the [submit.ncbi.nih.gov](http://submit.ncbi.nih.gov).
![NCBI submit Front](../img/ncbi_submit_genome_GWSS_1.png)

You'll need to know how your assembly was constructed and what assembly tools and dates or version of the tools. You will also need to have already created a BioProject and BioSample for this genome. You should also plan to deposit the raw sequence reads in the SRA but you can do that after the genome is deposited.

Depth of coverage of a genome calculated from total read amount and genome assembly size (bbcount or calculating this in mosdepth if you have a bam/cram file for your datasets aligned back to the genome).

Here's an example of calculating depth of coverage
```
module load mospdeth
mosdepth -x -n -t 16 remap_results Genome_remapped_reads.bam
```
and then run this R code to compute the mean depth of coverage. You can also add some additional filtering steps if you want to remove high coverage contigs from the average calculation, etc.

```R
library(dplyr)
mosd <- read.csv("remap_results.mosdepth.summary.txt",sep="\t",header=T)

mean(mosd$mean)
summary(mosd)
# remove high coverage for example using 100x as the cutoff in this example
nohighcov <- r %>% filter(mean < 100)
mean(nohighcov$mean)
summary(nohighcov)

```
This will give you a readout of the average depth of coverage across the contigs/chromosomes so you can fill in a number for depth of coverage in the genome submission entry.

Additional tools to calculate this include bbmap.sh
```bash
#SBATCH -p short -N 1 -n 48 --mem 64gb --out bbmap.log
module load BBMap
MEM=64
bbmap.sh -Xmx${MEM}g ref=GENOME.fasta in=LEFTREADS.fastq.gz in2=RIGHTREADS.fastq.gz statsfile=remapped_reads.bbmap_summary.txt
```

Once you have an annotations

### Uploading Genome data

The one, or many if you have a large genome the .sqn file will be split up by funannotate, sqn file(s) need to be uploaded to NCBI via either FTP or aspera. Aspera commands as shown before include a key that comes with aspera command line tool install as well as a private folder name for your account on NCBI which you can reveal in the Aspera command line instructions section of the page.

The upload step will look something like this on the UCR HPCC where *YOUREMAIL_SECRETFOLDER* will be what NCBI assigned to your account.

```
module load aspera
ascp -i $ASPERAKEY -QT -l 300m -k1 -d ./*.sqn subasp@upload.ncbi.nlm.nih.gov:uploads/YOUREMAIL_SECRETFOLDER/A_UNIQUE_NAME_FOLDER/
```
You specify a folder name (`A_UNIQUE_NAME_FOLDER`) for a particular project upload.



## Annotated genomes

For Bacteria there is a simple way to get these annotations and genome deposited in one go. Simply use the [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/) procedure.  The Genome submission portal is still used for uploading but then an option to annotate

If you would like to annotate a genome with PGAP before depositing you can also download and run [PGAP locally](https://github.com/ncbi/pgap).

# Assembled Transcriptomes

Use the Transcriptome Sequence Archive (TSA). These will still require a minimum length sequence filtering.

# Metagenome-assembled genomes

Instructions are [here](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/#metagen). First, create a BioProject and then create a Biosample for any host or environmental 'metagenome' samples of origin (the "physical" biosamples) as decribed above. Next email genomes@ncbi.nlm.nih.gov with a list of proposed taxonomic names for each MAG at its lowest rank for NCBI to review. 

> Here are the MAG names I sent to NCBI for consideration: 
> * Wolbachia sp. 
> * Candidatus Baumannia cicadellinicola  
> * Candidatus Sulcia muelleri 


> They sent back:
> * Wolbachia endosymbiont of Homalodisca vitripennis
> * Candidatus Baumannia cicadellinicola
> * Candidatus Sulcia muelleri

Use the NCBI provide names for MAGs to now create "organism"-specific biosamples. When creating the biosamples, choose "MIMAG Metagenome-assembled Genome". Include the BioProject ID made earlier during submission. Any Bin IDs can be included as the isolate name. Finally add a custom attribute 'derived-from' and add 'This BioSample is a metagenomic assembly obtained from the host XXXXXXXXXX BioSample: SAMNXXXXX.' to connect the MAG back to the physical biosample(s). If you've used a combined assembly (or co-assembly) approach, then instead here you will list all BioSamples such that 'This BioSample is a metagenomic assembly obtained from the metagenome BioSamples: SAMNXXX, SAMNXXX, and SAMNXXX.' You could also add custom attributes at this stage of % completion and % contamination for your MAGs. 

Finally, the MAGs can now be deposited via the WGS submision portal as genome assemblies (as described above) using the BioProject ID that we made earlier and the organism-specific BioSample IDs for each of the MAGs. For unannotated archaeal/bacterial MAGs, you only need a fasta file to submit these to NCBI. At the end of the submission process, you will be asked if you want to use NCBI's in-house annotation tool to annotate these, PGAP (as described above). Alternatively, you can provide .sqn files with your own annotation during deposition. 

Note - if you have a co-assembly you plan to deposit, in addition to "physical" and "organism" biosamples, you will also need to create a single BioSample that represents the 'combination of physical samples', SAMNXXX-SAMNXXX and include a note indicating that its a co-assembly such that 'This sample group is the combination of the XX individual BioSamples:AMNXXX, SAMNXXX, and SAMNXXX.' Then proceed as follows for a WGS submission.


