# Data Retrieval

## Downloading from NCBI Genomes

See how to use NCBI datasets to download
* https://github.com/1KFG/NCBI_fungi

## Downloading from FungiDB

* curl and wget examples from https://fungidb.org/fungidb/app/downloads/Current_Release/

## Downloading from JGI

### Using Globus

### Using scripts to parse download XML

(expand examples from https://github.com/zygolife/LCG_Assembly/tree/master/data/download_scripts
and https://github.com/stajichlab/DrylandAlgae_Genomes)


## Downloading from Ensembl Genomes

lftp and curl automation of download from https://ensemblgenomes.org/ see also
https://github.com/1KFG/2019_dataset/tree/master/scripts
and
https://github.com/1KFG/2019_dataset/

## Downloading from SRA

The SRAtoolkit provides built-in tools for downloading datasets.

This uses the HPCC module `workspace/scratch` which defines the `$SCRATCH` env variable.

```
#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -n 1 --mem 8gb --out dump.%a.log
module load sratoolkit
module load aspera
module load workspace/scratch
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
ACC=acc.txt
ACCNO=$(sed -n ${N}p $ACC | cut -f1)
ASPERA=$(which ascp)
if [ ! -f ${ACCNO}_1.fastq.gz ]; then
	prefetch --ascp-path "$ASPERA|$ASPERAKEY" $ACCNO
	fastq-dump --tmpdir $SCRATCH --defline-seq '@$sn[_$rn]/$ri' --split-files $ACCNO
	perl -i -p -e 's/_reverse//' ${ACCNO}_2.fastq
	perl -i -p -e 's/_forward//' ${ACCNO}_1.fastq
	pigz ${ACCNO}_[12].fastq

fi

```

### Parallelize fastq-dump

The SRA format to fastq conversion is applied by the fastq-dump tool.

Here's an example using array-jobs in slurm. There is a file called `lib/sra.txt` where each line is an SRA accession number.

This will use $SCRATCH variable (/scratch/$USERID/$SLURMJOBID) to write files so as long as there is enough space to run multiple downloads simultaneously and run the fastq dumping this can all proceed faster.  This example is intended to be an array job where each SRA download is run as a separate job so you would submit it as
`sbatch -a 1-10 download.sh` if the script was named `download.sh` and the `lib/sra.txt` had 10 lines for 10 SRR runs to download.

```
#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 16 --mem 8gb --out logs/download_sra.%a.log

module load parallel-fastq-dump
module load workspace/scratch

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi
SRAFILE=lib/sra.txt
FOLDER=input

MAX=$(wc -l $SRAFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SRAFILE"
  exit
fi
if [ ! -s $SRAFILE ]; then
	echo "No SRA file $SRAFILE"
	exit
fi
SRA=$(sed -n ${N}p $SRAFILE)
if [ ! -s ${SRA}_1.fastq.gz ]; then
	parallel-fastq-dump -T $SCRATCH -O $FOLDER  --threads $CPU --split-files --gzip --sra-id $SRA
fi
```
