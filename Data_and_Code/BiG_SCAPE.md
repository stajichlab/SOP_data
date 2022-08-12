#Using BiG-SCAPE to visualize secondary metabolites on a larger scale

I have created a template folder you can copy to the directory you want to work in. Check /shared/projects/BiG-SCAPE_template/.

A couple of things to prep before running this program:
* Make sure you have run the funannotate program, specifically antismash, on all your samples as you need the antismash gbk files to run this program.
* When running antismash make sure you have a locus number or a means to rename your GBK files to a different name, otherwise it will be labeled scaffold.
* Create a folder to store all your GBK files, I usually call it AntismashDir/ but it's up to you.
* Copy GBK files from your antismash folder into this folder. Make sure to only grab the smaller scaffold pieces and not the larger GBK file that is also found in that folder.
* Open the run_bigscape.sh file and see the basic command on the bottom to run BiG-SCAPE. Make sure the sure the name of the directory is the same as yours, if not change it.
* Run using sbatch and wait for results!


Good luck
-T
