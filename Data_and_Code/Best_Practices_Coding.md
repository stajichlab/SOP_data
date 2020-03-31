# Best Practices for Data management in Stajich lab

File naming

On HPCC the data are stored in
`/bigdata/stajichlab/` and within that folder is your own stuff. `bigdata`

Data organization
 - `bigdata` - your own
 - `shared/projects` are intended to be collaborative analysis environments.

 Sequence Data
 - primary sequence data is stored in `shared/SeqData` these are primary sequence from the sequencer and usually organized by date or flow cell code. They are to immutable and you analysis pipelines should link to this folder but not edit or update these files.

 - Best Practices would be to deposit the data at [NCBI SRA](NCBI_deposit) as soon as you know the data are accurate.


__more__
