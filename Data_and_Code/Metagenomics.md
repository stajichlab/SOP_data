# Stajich lab Metagenomics Best Practices
With all the new bioinformatics tools available it is hard to choose just one set to do your analysis with! 

Here is just a list (with scripts embedded) on the best practices held in our lab.
The series of steps include:
1. Trimming/removing index adapters
2. Normalizing read depth/coverage
3. Merging reads
4. Read Error Correction
5. Assembly MAGs
6. Map Reads to Assembly
7. Contig Binning Using AutoMeta
8. Check Assembly Quality (CheckM)
9. Confirm Identity GTDB-TK
10. 16S/ITS Reconstruction for Identity confirmation

After submitting samples to be sequenced you will receive demultiplexed samples. They should look like R1_001_clean1 or R2_001_clean1.

## 1. Trimming/removing index adapters
We need to remove the sequencing bits off our DNA.
```
#!/usr/bin/bash
#SBATCH -p batch --nodes=1 --ntasks=1 --cpus-per-task=12 --mem-per-cpu=24G --time=10-00:00:00 --job-name="MAG ASSEMBLY.2" --out logs/MAG.asm.1.%a.log

module load BBMap
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
	bbduk.sh -Xmx1g in1=${SAMPLE}_R1_001.fastq.gz in2=${SAMPLE}_R2_001.fastq.gz out1=${SAMPLE}_R1_CLEANEST.fastq out2=${SAMPLE}_R2_CLEANEST.fastq rcomp=t ktrim=r k=23 mink=11 ref=adapters hdist=1 maq=10 minlen=51 trimq=20 tpe tbo
done
```

### 1b. Confirm indexes are Removed - running FastQC
Run this on your original sample, and again after trimming to make sure your samples have been cleaned
```
module load fastqc

fastqc ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz --noextract -o FastQC

```

## 2. Normalizing Read Depth/Coverage
We need to correct for sequencing read depth.
```
#!/usr/bin/bash
#SBATCH -p batch --nodes=1 --ntasks=1 --cpus-per-task=12 --mem-per-cpu=24G --time=10-00:00:00 --job-name="MAG ASSEMBLY.2" --out logs/MAG.asm.2.%a.log

module load BBMap
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
	bbnorm.sh in1=${SAMPLE}_L001_R1_CLEANEST.fastq in2=${SAMPLE}_L001_R2_CLEANEST.fastq out1=${SAMPLE}_L001_R1_normalized.fastq out2=${SAMPLE}_L001_R2_normalized.fastq target=100 min=5
done    
```
## 3. Merging Your Reads
Let's put them all together now.
```
#!/usr/bin/bash
#SBATCH -p highmem --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=500G --time=10-00:00:00 --job-name="MAG ASSEMBLY.3" --out logs/MAG.asm.3.%a.log

module load BBMap
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
	bbmerge-auto.sh -Xmx500g in1=${SAMPLE}_L001_R1_normalized.fastq in2=${SAMPLE}_L001_R2_normalized.fastq out=${SAMPLE}_normalized_merged.fastq outu=${SAMPLE}_normalized_unmerged.fastq ihist=${SAMPLE}_ihist.txt ecct extend2=20
done
```

## 4. Read Error Correction
Correct our error reads
```
#!/usr/bin/bash
#SBATCH -p batch --nodes=1 --ntasks=1 --cpus-per-task=12 --mem-per-cpu=24G --time=10-00:00:00 --job-name="MAG ASSEMBLY.4" --out logs/MAG.asm.4.%a.log

module load spades/3.15.4
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
	spades.py -1 ${SAMPLE}_L001_R1_normalized.fastq -2 ${SAMPLE}_L001_R2_normalized.fastq --merged ${SAMPLE}_normalized_merged.fastq -o ${SAMPLE}_spades_error_corrected_reads -t 4 --meta --only-error-correction
done
```

## 5. Assemble MAGs using Fav Assembler - we use MetaSPAdes

```
#!/usr/bin/bash
#SBATCH -p batch --nodes=1 --ntasks=1 --cpus-per-task=12 --mem-per-cpu=24G --time=10-00:00:00 --job-name="MAG ASSEMBLY.5" --out logs/MAG.asm.5.%a.log

module load spades/3.15.4
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
	metaspades.py -1 ${SAMPLE}_L001_R1_normalized_error_corrected.fastq.gz -2 ${SAMPLE}_L001_R2_normalized_error_corrected.fastq.gz --merged ${SAMPLE}_normalized_merged_error_corrected.fastq.gz -s ${SAMPLE}_L001_R_unpaired_error_corrected.fastq.gz -o ${SAMPLE}_assembly_SPades --meta -t 4 --only-assembler
done
```

## 6. Map Reads to Assembly - Different tools than Autometa
```
#!/usr/bin/bash
#SBATCH -p batch --nodes=1 --ntasks=1 --cpus-per-task=12 --mem-per-cpu=24G --time=10-00:00:00 --job-name="MAG ASSEMBLY.6" --out logs/MAG.asm.6.%a.log

module load bwa-mem2/2.0
module load samtools/1.11
module load bwa/0.7.17
METAG=samples.csv  #this will be a separate file in your folder that lists your metagenome files if you have more than 1.

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

IFS=,
tail -n +2 $METAG | sed -n ${N}p | while read SAMPLE
do
# Build the index to map reads back
    bwa index ${SAMPLE}_contigs.fasta
# Map reads back using local alignment    
    bwa mem ${SAMPLE}_contigs.fasta ${SAMPLE}_L001_R1_normalized_error_corrected.fastq.gz ${SAMPLE}_L001_R2_normalized_error_corrected.fastq.gz -t 8 > ${SAMPLE}_contigs_aln_pe.sam
# View, converts SAM to BAM file, sort, index, and produce some coverage details
    samtools view -S -b ${SAMPLE}_contigs_aln_pe.sam > ${SAMPLE}_unsort_bwa.bam  
    samtools sort ${SAMPLE}_unsort_bwa.bam -o ${SAMPLE}_sorted.bam 
    samtools index  ${SAMPLE}_sorted.bam ## indexes BAM file
    samtools flagstat -@ 8 -O tsv ${SAMPLE}_sorted.bam > ${SAMPLE}_sorted_stats.tsv
    samtools coverage  ${SAMPLE}_sorted.bam -o ${SAMPLE}_sorted_coverage.tsv
    samtools depth  ${SAMPLE}_sorted.bam -o ${SAMPLE}_sorted_depth.tsv
done
```

Once assemblies are completed - you'll have final documents you can work from. Place this into a folder called genomes to progress.
## 7. Autometa Steps
Using Autometa to help take our assembly reads, determine taxonomy, predictions, etc.

cat 01_make_cov.sh
```
#!/usr/bin/bash
#SBATCH -p batch --mem 32gb -N 1 -n 24 --out logs/autometa_taxonomy.%a.%A.log

# see module load below

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi

DATABASES=databases
if [ ! -d $DATABASES ]; then
    ln -s /srv/projects/db/autometa/1.0.2 $DATABASES
fi
DAT=samples.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SAMPLE
do
    echo $SAMPLES
    GENOMEFILE=genomes/${SAMPLE}.dna.fasta
    COVTAB=coverage/$SAMPLE.coverage.tab
    AUTOMETAOUT=autometa_runs
    mkdir -p $AUTOMETAOUT
    if [ ! -f $COVTAB ]; then
	bash pipeline/autometa/00_make_cov.sh $N
    fi
    if [ ! -d $AUTOMETAOUT/$SAMPLE ]; then
	module load autometa/1.0.2
	make_taxonomy_table.py -a $GENOMEFILE -p $CPU -o $AUTOMETAOUT/$SAMPLE --cov_table $COVTAB
    fi
done
```
cat 02_combine_contig_cov_table.sh

```
#!/usr/bin/bash
#SBATCH -N 1 -n 1 -p short --out logs/contig_cov_table.%a.log

module load autometa/1.0.2
source activate autometa
module load git

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	echo "need to provide a number by --array or cmdline"
	exit
    fi
fi

DAT=samples.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SAMPLE
do
    echo $SAMPLE
    OUT=coverage
    GENOMEFILE=genomes/${SAMPLE}.dna.fasta
    COVTAB=$OUT/$SAMPLE.coverage.tab
    CONTIGCOV=$OUT/$SAMPLE.contig_cov_table.tab
    if [ ! -f $CONTIGCOV ]; then
	make_contig_table.py -a $GENOMEFILE -c $COVTAB -o $CONTIGCOV
    fi
done
```
cat 03_run_autometa.sh

```
#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/autometa_run.%a.%A.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	echo "need to provide a number by --array or cmdline"
	exit
    fi
fi
AUTOMETAOUT=autometa_runs
COVERAGE=coverage
DAT=samples.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SAMPLE
do
    echo "Processing $SAMPLE"
    COVTAB=$COVERAGE/$SAMPLE.coverage.tab

    if [ -d $AUTOMETAOUT/$SAMPLE ]; then
	module load autometa/1.0.2
	module load git
	time run_autometa.py -k bacteria -a $AUTOMETAOUT/$SAMPLE/Bacteria.fasta --ML_recruitment \
		--processors $CPU --length_cutoff 1500 --taxonomy_table $AUTOMETAOUT/$SAMPLE/taxonomy.tab -o $AUTOMETAOUT/$SAMPLE/Bacteria_run -v $COVTAB
    fi
done
```

## 8. Check Assembly Quality (CheckM)
We have assembled and built our metagenomes. Now how do we know we have a full genome and not a small subsample? We run CheckM.

```
#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/checkm_run.%a.%A.log

module load checkm
module unload miniconda2/4.4.10
module load anaconda3/4.5.4

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

checkm taxonomy_wf phylum Cyanobacteria -t 8 -x faa --genes annotated_protein_genomes/ checkM_run/checkM/ -f checkM_taxonomy_wf_results --tab_table
```
Notice how the above script has it listed for Cyanobacteria, it depends on your samples, do you know what type of bacteria has come out? (Check your taxonomy tables to confirm for each strain). Find in checkm documentation what would your phylum be and replace in script to run.

## 9. Confirm Identity GTDB-TK
You have your compiled MAGs but want to confirm which species your organisms are. Running GTDB-TK will place your strains in a giant phylogenetic tree. This can be computationally taxing to visualize. See if you have the specifications to run this tree vizualizer (we use FigTree). 

cat 01_identify_gtdbtk.sh 

```
#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/gtdbtk_run_identify.%a.%A.log


CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

module load gtdbtk

gtdbtk identify --genome_dir genomes/ --out_dir GTDB-TK/identify/ --extension fasta --cpus 2
```

02_align_gtdbtk.sh

```
#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/gtdbtk_run_align.%a.%A.log


CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

module load gtdbtk

gtdbtk align --identify_dir GTDB-TK/identify/ --out_dir GTDB-TK/align/  --cpus 2
```

cat 03_classify_gtdbtk.sh
```
#!/usr/bin/bash
#SBATCH -p short --mem 200G  -n 24 --out logs/gtdbtk_classify_run.%a.%A.log


CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

gtdbtk classify --genome_dir GTDB-TK/genomes/ --align_dir GTDB-TK/align/ --out_dir GTDB-TK/classify/ --extension fasta --cpus 2
```

## 10. 16S/ITS Reconstruction for Identity confirmation
Another method of pulling out the identity of your MAGs is to use PhyloFlash
In-progress ...


# SnakeMake Pipeline in Progress ...
## Install the snakemake pipeline from conda package [**metagenome-atlas**](https://metagenome-atlas.readthedocs.io/en/latest/usage/getting_started.html#setup)
1. First, create an environment which use **python3.8**.

    `conda create -n metagenome-atlas python=3.8`

    Note than some of the dependencies is not **python3.9** supported so be sure to use **python3.8**. After create the environment, activate it by:
    
    `conda activate metagenome-atlas`

2. Install the package **mamba**, which is a faster version of **conda**. 

    `conda install -c conda-forge mamba`
    
    After **mamba** being installed, you can later switch from `conda install [package]` to `mamba install [package]`. The speed is way faster than conda.

3. Next, install the package **metagenome-atlas** through **mamba**.
    
    `mamba install metagenome-atlas`
    
    It is recommended to use **mamba** to install **metagenome-atlas**. Although **metagenome-atlas** can be found on [anaconda](https://anaconda.org/bioconda/metagenome-atlas), you might have chance to run into dependency troubles while install it directly through conda.
    
4. Then you can execute the `atlas --help` to show the metagenome-atlas helping page
    
    ```
    Usage: atlas [OPTIONS] COMMAND [ARGS]...

      ATLAS - workflows for assembly, annotation, and genomic binning of
      metagenomic and metatranscriptomic data.

      For updates and reporting issues, see: https://github.com/metagenome-
      atlas/atlas

    Options:
      --version   Show the version and exit.
      -h, --help  Show this message and exit.

    Commands:
      download  download reference files (need ~50GB)
      init      prepare configuration file and sample table for atlas run
      run       run atlas main workflow
    ```
