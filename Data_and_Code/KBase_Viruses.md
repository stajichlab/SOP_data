# Using KBase for viral discovery and annotation from metagenomes

KBase is an online web tool from the Department of Energy that is focused on providing bioinformatics infrastructure for genomics and metagenomics related anlayses. We have found that this tool works well for undergraduate researchers and those new to bioinformatic analyses as a stepping stone to command line work.


KBase can be accessed by creating an account or signing in through one's ORCID. Within KBase the [Stajich Lab](https://narrative.kbase.us/#/orgs/stajichlab) has a group that one can share and associate narratives with. Narratives are essentially IPython notebooks (or bioinformatic labratory notebooks) which contain all of the KBase modules that are run and what data is used as well as any resulting outputs. KBase narratives can also be 'forked' or copied and then edited allowing for users to re-use or edit a narrative as a jumping off point for their own analyses.


The viral protocol that we follow is currently provided as a [KBase tutorial narrative 'Viral End-To-End Analysis'](https://kbase.us/n/75811/85/) which is the current protcol we use in lab. Dr. Benjamin Bolduc also gives a webinar on use of this narrative through KBase which is on [YouTube](https://www.youtube.com/watch?v=WbQeOzdyTbc) and is highly recommended.


An almost identical protcol is on [Protocols.io](https://www.protocols.io/view/processing-a-viral-metagenome-using-ivirus-buisnuee), but which uses CyVerse which our undergraduates found harder to work with.


Notes on the protcol: 
1. VirSorter will sort predicted viruses into 6 categories (1-3 for predicted lytic/viral, and 4-6 are prophage that are integrated). Categories 1-2 and 4-5 are the most confident are most people seem to not use the remaining categories (3/6). 
2. VContact2 does not automatically assign taxonomy to viruses - the user must look at the resulting clustering of any viral data provided with the reference database chosen and then assess based on these clusters (and follow-up work) whether or not to assign taxonomy to each virus. 
3. Currently, there is no way to combine assembly objects in KBase and VContact2 will not take in Assembly Set objects - thus - VContact2 can only take in one assembly at a time and any resluting comparative work wil need to be done off of KBase (e.g. constructing vOTUs based on clustering etc to look at viral abundance ACROSS samples). The only current work around is to combine read data at the beginning of the pipeline into one set and do a co-assembly - though depending on the data size, you might run into memory errors here.
4. Currently you need to use the VContactv2 beta application if you want to download the resulting network for visualization in Cytoscape (see Protocols.io tutorials for Cytoscape instructions). 
