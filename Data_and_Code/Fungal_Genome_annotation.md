# Fungal Genome Annotation

In our lab we typically use [Funannotate](https://funannotate.readthedocs.io/) for annotation though [MAKER](https://yandell-lab.org/software/maker.html) is also another pipline that we support and have used in the past.

_Additional tutorial material will be provided here in an update_

For now see:
https://github.com/stajichlab/funannotate_template


# Running Gene prediction Training

Using RNA-seq data Funannotate can improve the gene prediction parameters for SNAP and Augustus.

If you are utilizing RNASeq for training and improvement of gene models you can take advantage of the extra speed up that running PASA with a mysql server (instead of the default, SQLite).  To do that on HPCC this requires seting up a mysql server using singularity package.

## Setting up MariaDB/Mysql

You need to create a file `$HOME/pasa.CONFIG.template` this will be customized for your user account. Copy it from the system installed PASA folder.
A current version on the system is located in `/opt/linux/centos/7.x/x86_64/pkgs/PASA/2.3.3/pasa_conf/pasa.CONFIG.template`.
Doing `rsync /opt/linux/centos/7.x/x86_64/pkgs/PASA/2.4.1/pasa_conf/pasa.CONFIG.template ~/`
This can also be done automatically with the latest version of PASA on the system
```bash
module load PASA/2.4.1
FOLDER=$(dirname `which pasa`)
rsync -v $FOLDER/../pasa_conf/pasa.CONFIG.template
```

You will need to edit this file which has this at the top. The MYSLQSERVER part will get updated by the mysql setup step later so leave  it alone.
You will need to fill in the content for `MYSQL_RW_USER` and `MYSQL_RW_PASSWORD` too.

```
# server actively running MySQL
# MYSQLSERVER=server.com
MYSQLSERVER=localhost
# Pass socket connections through Perl DBI syntax e.g. MYSQLSERVER=mysql_socket=/tmp/mysql.sock

# read-write username and password
MYSQL_RW_USER=xxxxxx
MYSQL_RW_PASSWORD=xxxxxx
```

On the UCR HPCC [here are directions](https://github.com/ucr-hpcc/hpcc_slurm_examples/tree/master/singularity/mariadb) on how to setup your own mysql instance in your account using [singularity](https://sylabs.io/docs/). If you were running funannotate on your own linux/mac setup you would just do a native mysql/mariadb install and have the server running on your local machine.

The HPCC instructions include the steps to initialize a database followed by you will start a job that will be running which has the mysql instance. This db server will need to be started before you start annotating and be shutdown when you are finished. I often give it a long life like 2 weeks but it can be stopped at any point too.

## If you are to run this in our lab setting (stajichlab) on the cluster follow the directions below.
Make sure your directory has the following:
* genomes/    This folder should contain the genomes you would like to annotate. Generally this is a genome you have run AAFTF or any other assembly program on. We recommend you also sort this genome from largest scaffold to the smallest too. You can do this using [AAFTF sort](https://github.com/stajichlab/AAFTF).
* logs/    This will be a folder where all the log files will be added to. This way you can come to this folder and see what went wrong!
* pipeline/    This folder will hold your scripts. You can find [them here](https://github.com/stajichlab/funannotate_template/pipeline)
* samples.csv    This text file will contain all your strain information. The [scripts here](https://github.com/stajichlab/funannotate_template) are prepared for more than one genome to be annotated in an array. If you don't have more than one sample you can still add this into the text file.
* lib/augustus/3.3/config/    Folder you can find this.. (WHERE??) This folder can also contain your protein.aa file for your genome, along with your [.sbt file](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) you should generate from NCBI if you are planning on submitting this genome online.
Your samples.csv file should look like:
```
SPECIES,STRAIN,PHYLUM,LOCUS
Fusarium_odoratissimum,CF115,Sordariomycetes,CF115
```
All other folders will be created by running the steps in the pipeline folder.

The Funannotate steps look like:
1. Run MySQL database. (000_start_mysql.sh)
2. Build a Repeat Modeler database on your genome(s). (00_repeatmodeler.sh)
3. Run Repeat Masker on your genome(s). (01_mask_denovo.sh)
4. IF you have RNA-Seq data, you can run the training step, if not you can skip this step. (02_train_RNASeq.sh)
5. Run the prediction step using Augustus. (03_predict.sh)
6. Run the update step (if you have run the training step) to add the 5' and 3' UTRs. (04_update.sh)
7. Run the antismash step to predict secondary metabolites in your genome(s). (05a_antismash_local.sh)
8. Run the iprscan step to predict functional protein and GO terms for your predicted proteins. (05b_iprscan.sh)
9. Run the final Annotation step in Funannotate that adds all the information generated from above into final versions of protein, sbt, annotations.txt, gbk and gff3 files. (06_annotate_function.sh)

Paste the following MySQL code into your pipeline directory. Submit this to sbatch then commence through the rest of the Funannotate steps. There is no need to change anything in this script.

```
cat 000_start_mysql.sh
```
[000_start_mysql.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/000_start_mysql.sh)
This script will run for 5 days, make sure you finish all the annotation steps in this time.

## Next pre-step is to run RepeatModeler.
This program builds a repeated element database out of your genome in order to be masked properly.

```
cat 00_repeatmodeler.sh
```

```
#!/usr/bin/bash -l
#SBATCH -p batch --time 2-0:00:00 --ntasks 8 --nodes 1 --mem 120G --out logs/repeatmodeler_attempt.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=genomes

mkdir -p repeat_library

SAMPFILE=samples.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $(expr $MAX) ]; then
    MAXSMALL=$(expr $MAX)
    echo "$N is too big, only $MAXSMALL lines in $SAMPFILE"
    exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUS
do
  name=$(echo -n ${SPECIES}_${STRAIN} | perl -p -e 's/\s+/_/g')
  echo "$name"
     module unload perl
     module unload python
     module unload miniconda2
     module unload anaconda3
     module load RepeatModeler
     module load ncbi-blast/2.13.0+
     export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
	#makeblastdb -in $INDIR/$name.sorted.fasta -dbtype nucl -out repeat_library/$name
	BuildDatabase -name repeat_library/$name $INDIR/$name.sorted.fasta
	RepeatModeler -database repeat_library/$name -pa $CPU
	#-LTRStruct
done
```

In order to run this step properly you must submit like:
```
sbatch --array=1 pipeline/00_repeatmodeler.sh
```
If you have more than one genome, make sure your samples.csv file has this information inside, and the genome can be found in the genomes folder and the way to run this looks like:
```
sbatch --array=1-5 pipeline/00_repeatmodeler.sh
```

## Once you have created a repeat modeler database, you can now mask using repeatmasker.

```
cat 01_mask_denovo.sh
```
[01_mask_denovo.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/01_mask_denovo.sh)

You would be submitting this script just like above:
```
sbatch --array=1-5 pipeline/01_mask_denovo.sh
```
If this was run successfully, your genomes/ folder will have your original genome along with a masked.genome found here.
Make sure your genome ends with GENUS_SPECIES.sorted.fasta in order to run this properly.

## Once your genome is masked properly you can now run the training step. This step is not required.
If you have RNA-seq data on your genome, you can run this step. Otherwise you should skip to the next step (03_predict) in the Funannotate pipeline.

```
cat 02_train_RNASeq.sh
```
[02_train_RNASeq.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/02_train_RNASeq.sh)

This portion of the pipeline requires a couple more items to be run properly. Make sure you have the lib/RNASeq folder here. If you have RNA-Seq data available for your genome, or even RNASeq data of a close species of your organism you would like to use to better predict your genome, you should add them here.
Your folder should look like.
lib/RNASeq/YOUR_GENOME_SPECIES_NAME/
Inside this folder should have two files Forward.fq.gz and Reverse.fq.gz. These two files you can find through SRA and download into this folder and rename to these two names.
Once complete you should see a annotate/YOUR_GENOME_SPECIES/training folder.

## Once completing the RNA training step, you can run the predict step.
```
cat 03_predict.sh
```
[03_predict.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/03_predict.sh)
```
Now that you have finally set up and masked and trained your genome, running this step is one of the main results from this pipeline. You can either use the RNASeq data to generate a better predicted genome, or you can skip the prior step and move on to this step to run the prediction denovo.
Make sure the line in this script says: ```tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG```
Check your annotate/YOUR_GENOME_SPECIES/predict_results to see if the program has run to completion. You will see a series of files including gbk, gff3 and protein files. But you are not done just yet!!

## If you have run the RNASeq step then you can run the next step which is the update step in the Funannotate program.
```
cat 04_update.sh
```
[04_update.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/04_update.sh)
Make sure the line in this script says: ```tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG```
Check your annotate/YOUR_GENOME_SPECIES/update_results to see if the program has run to completion. You will see a series of files including gbk, gff3 and protein files. But you are not done just yet!!

## Now you can run two alternative annotation programs, [Antismash](https://fungismash.secondarymetabolites.org/#!/start) and [IPRScan](https://www.ebi.ac.uk/interpro/search/sequence/).

```
cat 05a_antismash_local.sh
```
[05a_antismash_local.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/05a_antismash_local.sh)
Make sure the line in this script says: ```tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG```
```
cat 05b_iprscan.sh
```
[05b_iprscan.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/05b_iprscan.sh)
Make sure the line in this script says: ```tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG```

At the end of these two scripts you should have folders.
* antismash_local
* annotate_misc

Make sure you look through these folders and see if annotate/YOUR_GENOME_SPECIES/antismash_local/ has gbk files, there should be one for each scaffold. While the annotate/YOUR_GENOME_SPECIES/annotate_misc folder should have an iprscan.xml file (and it should NOT be empty).

# FINAL Step for Funannotate is to run the annotation step.
```
cat 06_annotate_function.sh
```
[06_annotate_function.sh](https://github.com/stajichlab/funannotate_template/blob/main/pipeline/06_annotate_function.sh)
Make sure the line in this script says: ```tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG```
Once the script has run successfully you should find annotate/YOUR_GENOME_SPECIES/annotate_results and within you will find gbk, gff3, protein and a couple of other file types once complete.

Congratulations! You have run the Funannotate pipeline!

-T
