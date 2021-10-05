# Data Retrieval

## Downloading from NCBI Genomes

## Downloading from FungiDB

## Downloading from JGI

## Downloading from Ensembl Genomes

## Downloading from SRA

The SRAtoolkit provides built-in tools for downloading datasets.

This uses the HPCC module `workspace/scratch` which defines the `$SCRATCH` env variable. 

```
#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 1 --mem 8gb --out dump.%a.log
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

```
#!/usr/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 8gb --out logs/download_sra.%a.log

module load sratoolkit
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
