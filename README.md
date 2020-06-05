# VEP variant annotation pipeline

## 1. Installation
This section describes how to set up VEP, download all necessary cache files, and install LoFtee plugin.

### 1.1. Setting up VEP

1. Load `Singularity` module:
   ```
   module load singularity
   ```

2. Build `Singularity` image:
   ```
   singularity build vep.sif docker://dtaliun/vep_samtools:latest
   ```
   This step may take around 1h.

3. Download VEP cache files into local `vep_cache` directory:
   ```
   mkdir `pwd`/vep_cache
   export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
   singularity run -B `pwd`/vep_cache:/opt/vep/.vep vep.sif INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep
   ```
   This step may take more than 1h.

### 1.2. Setting up LoFtee

More detailed instructions on how to set up LoFtee are [here](https://github.com/konradjk/loftee).

1. You must clone LoFtee repository into your local `vep_cache` directory:
   ```
   cd vep_cache
   git clone https://github.com/konradjk/loftee.git
   ```

2. Switch to the GRCh38 branch (or the branch you need):
   ```
   git checkout grch38
   ```
   
3. Download all necessary databases (based on human genome build you plan to use) as described [here](https://github.com/konradjk/loftee). BUT: all these files must be stored under your `vep_cache` directory e.g. `vep_cache/loftee_b38`.

### 1.3. Conclusion

After above steps, your local `vep_cache` directory should be similar to this:
```
|- vep_cache
   |- homo_sapiens (directory with VEP databases)
   |- loftee (loftee scripts)
   |- loftee_b38 (loftee databases)
```

## 2. Running

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/vep_pipeline.git
   ```

2. Modify `nextflow.config` configuration file.
     * `params.vcfs` -- path to your VCF/BCF file(s). You can use `glob` expressions to selecect multiple files.
     * `params.vep_cache` -- full path to your local `vep_cache` directory.
     * `params.vep_flags` -- flags you want to pass to VEP.
     * `params.loftee_dir` -- name of the directory within `vep_cache` with LoFtee scripts.
     * `params.loftee_db_dir` -- name of the directory within `vep_cache` with LoFtee database files.
     * `process.container` -- full path to the `Singularity` image file (see step 1.1.).
     * `executor.$slurm.queueSize` -- maximal number of SLURM jobs to submit at once.
  
3. Run pipeline:
   ```
   module load nextflow
   module load singularity
   nextflow run Annotate.nf -w ~/scratch/work_directory
   ```
   Important: when working on Compute Canada HPC, set working directory to ~/scratch/\<new directory name\>. This will speed up IO and also save space on your `project` partition. After the execution, if there were no errors and you are happy with the results, you can remove this working directory.
  
## 3. Known pitfalls

1. You may not be able to execute `nextflow` directly from the Compute Canada login nodes due to the 8Gb memory limit per user. One alternative is to start interactive slurm job and submit all commands from it e.g.:
   ```
   salloc --time=2:00:00 --ntasks=1 --mem-per-cpu=16G
   ```
   Or submit a batch job with nextflow command e.g.:
   ```
   module load nextflow
   module load singularity
   sbatch --time=2:00:00 --ntasks=1 --mem-per-cpu=16G --wrap="nextflow run Annotate.nf -w ~/scratch/work_directory"
   ```
   Make sure you specify enough time. VEP annotation is typically fast, but total `nextflow` execution time will depend on how busy the SLURM queue is.

2. Sometimes `nextflow` will crash with error `Failed to submit process to grid scheduler for execution`. Most probably the SLURM queue was too busy and thus slow to respond. Your results were not lost, just resume `nextflow` execution with the following command and `nextflow` will continue from where it finished:
   ```
   nextflow run Annotate.nf -w ~/scratch/work_directory -resume
   ```
