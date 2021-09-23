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