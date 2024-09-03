# Phylogenomics with PHYling

PHYling is a package for building species trees by identifying homologs to a core set of phylogenetic marker genes and performs the steps to identify, extract sequences, align, and build trees. It is mostly aimed around using existing protein profile Hidden Markov Models. Taking advantage of the work from OrthoDB and BUSCO team we can use markers from those projects as input.

## Fungi work

A set of protein coding markers were developed over the course of the **A**ssembling the **F**ungal **T**ree **O**f **L**ife (AFTOL). (Citations) The first set of 3 protein coding markers were RPB1, RPB2, and EF1alpha we call AFTOL1. Further work using a collection of **N** genomes identified a set of 70 markers called AFTOL2 which were all named by their _Saccharomyces cerevisiae_ gene names ([Floudas 2012](http://dx.doi.org/10.1126/science.1221748)). Further work from the Joint Genome Institute orthologous gene clustering of Dikarya fungal genomes in ~2011 generated a set of markers we named JGI_1086. This set of 434 were generally single copy genes found in the Dikarya fungi.

Additional work using the BUSCO / OrthoDB v9 release - fungi_odb9 - have also been deployed as markers for phylogenomic studies (cite Li et al) which number ~300 (Check that). Some of my analyses found these markers did not perform as well as the JGI_1086 set. The version 10 OrthoDB markers were released in (check date) and the fungi_odb10 set tested. Tests with this marker set on the 1000 fungal genome dataset has demonstrated better congruence with previous expectations of the fungal phylogenetic tree.

# Preparing data for PHYling

Checking out the code

Creating folders:
 - pep
 - cds

 config.txt file

# Running PHYling

PHYling can be downloaded from [github](https://github.com/stajichlab/PHYling_unified). This can be checked out with.
```bash
git clone https://github.com/stajichlab/PHYling_unified.git
```
### Tool Dependencies
 - [HMMer](http://hmmer.org/) 3.x including the easel tools (esl-XXX).
 - [TrimAl](http://trimal.cgenomics.org/)
 - [muscle](http://drive5.com) if you want to do _de novo_ alignment instead of aligning to the HMMs with hmmalign

### Data Dependencies
 The folder of HMM models need to be a directory called HMM. Each model set (eg JGI_1086, fungi_odb10) are in a folder with that name within the HMM folder.
 Several datasets are available from [1KFG HMMs repository](https://github.com/1KFG/PHYling_HMMs_fungi).
 - AFTOL_1
 - AFTOL_70
 - Roz200
 - JGI_1086

Additional current best practices markers are available from [BUSCO Fungi](https://busco.ezlab.org/frames/fungi.htm) [fungi_odb10] or for other organisms [BUSCO v4 / ODB10](https://busco.ezlab.org/busco_v4_data.html)

## Preparing input data

Protein files are typically named with a common suffix (`aa.fasta`) and are located in the `pep` folder. This default can be changed in the `config.txt` file with the `PEPDIR` and `PEPEXT` parameters.

Here's a complete `config.txt` example
```text
#basic setup needs to be changed
PHYLING_DIR=__FIXME__
HMM=__CHANGEME__
PREFIX=__CHANGEME__
HMM_FOLDER=HMM
PEPDIR=pep
CDSDIR=cds
INPEPEXT=aa.fasta
INCDSEXT=cds.fasta
ALLSEQNAME=allseq
OUTPEPEXT=aa.fa
OUTCDSEXT=cds.fa
LISTFILE=pepfile.lst
BESTHITEXT=best
HMMSEARCH_CUTOFF=1e-30
HMMSEARCH_OUTDIR=search
ALN_OUTDIR=aln
LANGUAGE=en
#job runs
JOBCPU=2 # per job CPUs
TOTALCPU=8 # total CPUs to use
QUEUEING=parallel
QUEUE=
TEMP=/tmp/__CHANGEME__
# tree building - not yet integrated into PHYling
OUTGROUP=__CHANGEME__
EXTRARAXML=
EXTRAIQTREE="-nt AUTO -m TESTMERGE -bb 1000 -alrt 1000"
```
## Initialize

```
PHYling init
```

## Searching for homologs

```bash
PHYling search
```

```bash
PHYling search -q slurm
```

```bash
PHYling search -q slurm --force
```

## Aligning  

```bash
PHYling aln
```

```bash
PHYling aln -q slurm
```
