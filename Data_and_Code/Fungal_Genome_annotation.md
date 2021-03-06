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
