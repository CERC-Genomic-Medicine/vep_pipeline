# VEP variant annotation pipeline

## Prerequisites

The following software is required:
- Apptainer (tested with version 1.2.4)
- Nextflow (tested with version 23)
-  **Only when using `enable_summary = true`**. Python 3 (tested with version 3.7.7) with the following packages: pysam, nbconvert, ipykernel, pandas.

## Workflow
> [!CAUTION]
> - Your input VCFs must be indexed and have corresponding `.tbi` files
> - Input VCFs from the same study must have the same prefix, which can't be a chromosome name.
> - Input VCFs can be split by chromosome as long as they all have the same prefix e.g., [prefix].chr1, [prefix].chr2, ..., [prefix].chrY.

![Execution diagram](Diagram.png)

## 1. Installation
> [!NOTE]
> If you are a member of the CERC-Genomic-Medicine team, then all necessary files to run VEP are already installed on the server. Consult with your colleagues on where to find them.

This section describes how to set up VEP, download all necessary cache files, and install LoFtee plugin.

### 1.1. Setting up VEP

1. Load `apptainer` module:
   ```
   module load apptainer
   ```

2. Build the SIF image with additional tools: `samtools`, `bcftools`, and `DBD::SQLite`.

   Use [SylabsCloud](https://cloud.sylabs.io/) free Remote Builder service to create the SIF image remotely.
   Specify the following in your definition file (i.e. `.def` file):
   ```
   Bootstrap: docker
   From: ensemblorg/ensembl-vep:latest
   
   %post
        apt-get update -y
        apt-get install -y samtools
        apt-get install -y bcftools
        apt-get install -y libdbd-sqlite3-perl
   ```

   After you built the SIF image using the Remote Builder web interface, pull the image to the cluster.

   Confirm that the `SylabsCloud` remote service is available (if not, you will need to add it):
   ```
   apptainer remote list
   ```

   Login to the SylabsCloud remote service
   ```
   apptainer remote login SylabsCloud
   ```
   
   Pull the SIF image using the repository name you speficied when creating the SIF image:
   ```
   apptainer pull vep.sif library://dtaliun/remote-builds/vep:23may2024
   ```

4. Download VEP cache files into local `vep_cache` directory:
   ```
   mkdir `pwd`/vep_cache
   export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
   apptainer run -B `pwd`/vep_cache:/opt/vep/.vep vep.sif INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep
   ```
   This step may take more than 1h.

### 1.2. Setting up LoFtee

More detailed instructions on how to set up LoFtee are [here](https://github.com/konradjk/loftee).

1. You must clone LoFtee repository into your local `vep_cache` directory:
   ```
   cd vep_cache
   git clone https://github.com/konradjk/loftee.git loftee_GRCh37
   git clone https://github.com/konradjk/loftee.git loftee_GRCh38
   cd loftee_GRCh38
   git checkout grch38
   cd ..
   ```
   
2. Download all necessary databases (based on human genome build you plan to use) as described [here](https://github.com/konradjk/loftee) into your `vep_cache` directory into folders `loftee_db_GRCh37` and `loftee_db_GRCh38`. These should include: GERP conservation scores (only for GRCh38), human_ancestor.fa files, SQL databases with PhyloCSF metrics (SQL files must be unzipped).

### 1.3. Setting up CADD

1. Download VEP plugins into `vep_cache` directory:
   ```
   cd vep_cache
   git clone https://github.com/Ensembl/VEP_plugins.git Plugins
   cd ..
   ```
2. Download CADD scores for GRCh37 and GRCh38 builds
   ```
   cd vep_cache
   
   mkdir CADD_GRCh37
   cd CADD_GRCh37
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz.tbi
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz.tbi
   cd ..
   
   mkdir CADD_GRCh38
   cd CADD_GRCh38
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi
   cd ..
   
   cd ..
   ```

### 1.4. Adding custom VEP plugins

1. Copy `Plugins/CONTEXT.pm` file from this repository to the `vep_cache/Plugins` directory.
   ```
   cd Plungins
   wget https://raw.githubusercontent.com/CERC-Genomic-Medicine/vep_pipeline/master/Plugins/CONTEXT.pm
   cd ..
   ```

### 1.5. Conclusion

After above steps, your local `vep_cache` directory should be similar to this:
```
|- vep_cache
   |- homo_sapiens (directory with VEP databases)
   |- loftee_GRCh37 (loftee scripts for build GRCh37)
   |- loftee_GRCh38 (loftee scripts for build GRCh38)
   |- loftee_db_GRCh37 (loftee databases for build GRCh37)
   |- loftee_db_GRCh38 (loftee databases for build GRCh38)
   |- Plugins (CADD plugin)
   |- CADD_GRCh37 (CADD scores for build GRCh37)
   |- CADD_GRCh38 (CADD scores for build GRCh38)
```

## 2. Running

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/vep_pipeline.git
   ```

2. Modify `nextflow.config` configuration file.
     * `params.vcfs` -- path to your VCF/BCF file(s). You can use `glob` expressions to selecect multiple files.
     * `params.assembly` -- set to "GRCh37" or "GRCh38".
     * `params.vep_cache` -- full path to your local `vep_cache` directory.
     * `params.vep_flags` -- flags you want to pass to VEP.
     * `params.loftee_flags` -- comma-separated list of additional LoFtee flags (with leading comma). Flags `loftee_path`, `gerp_bigwig`, `human_ancestor_fa`, and `conservation_file` are set automatically based on the selected `assembly`.
     * `enable_summary` -- set to `true` if you want to generate HTML summary files.
     * `process.container` -- full path to the `Singularity` image file (see step 1.1.).
     * `executor.$slurm.queueSize` -- maximal number of SLURM jobs to submit at once.
  
3. Run pipeline:
   ```
   module load nextflow
   module load apptainer
   nextflow run Annotation.nf -w ~/scratch/work_directory
   ```
   Important: when working on Compute Canada HPC, set working directory to ~/scratch/\<new directory name\>. This will speed up IO and also save space on your `project` partition. After the execution, if there were no errors and you are happy with the results, you can remove this working directory.
  
## 3. Custom VCFs

In this section, we will explain how to integrate custom VCF files, such as gnomAD v2, into your VEP command line using this pipeline.

1. Download the custom VCF into the cache. As an exemple, we downloaded the gnomAD v2 VCF file which is located at:
 ```
   /path/to/vep_cache/custom_vcf/gnomad.exomes.r2.1.1.sites.liftover_grch38.PASS.noVEP.vcf.gz
 ```

2. To integrate the custom VCF into your VEP command, add the relevant flags to the `nextflow.config` file under the `--custom` flag. Below is an example configuration:
```
vep_flags = "--sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --af_gnomade --af_gnomadg --pubmed --shift_hgvs 0 --allele_number --buffer_size 10000 --custom /path/to/vep_cache/custom_vcf/gnomad.exomes.r2.1.1.sites.liftover_grch38.PASS.noVEP.vcf.gz,gnomad_exomes,vcf,exact,0,AN_nfe,AC_nfe,non_cancer_AN_nfe,non_cancer_AC_nfe"
```
In this case, the custom flag integrates the gnomAD v2 VCF with the following parameters:

    path: Path to the custom VCF file.
    identifier: A name for the custom annotation (e.g., gnomad_exomes).
    type: The file type (e.g., vcf).
    match_type: How to match the data (e.g., exact).
    cols: Columns to include from the VCF (e.g., 0,AN_nfe,AC_nfe,non_cancer_AN_nfe,non_cancer_AC_nfe).


3. Run pipeline.
     Once you have updated your nextflow.config, run your VEP command as usual. The custom VCF annotations will be included in the output.

## 4. Known pitfalls

1. You may not be able to execute `nextflow` directly from the Compute Canada login nodes due to the 8Gb memory limit per user. One alternative is to start interactive slurm job and submit all commands from it e.g.:
   ```
   salloc --time=2:00:00 --ntasks=1 --mem-per-cpu=16G
   ```
   Or submit a batch job with nextflow command e.g.:
   ```
   module load nextflow
   module load apptainer
   sbatch --time=2:00:00 --ntasks=1 --mem-per-cpu=16G --wrap="nextflow run Annotation.nf -w ~/scratch/work_directory"
   ```
   Make sure you specify enough time. VEP annotation is typically fast, but total `nextflow` execution time will depend on how busy the SLURM queue is.

2. Sometimes `nextflow` will crash with error `Failed to submit process to grid scheduler for execution`. Most probably the SLURM queue was too busy and thus slow to respond. Your results were not lost, just resume `nextflow` execution with the following command and `nextflow` will continue from where it finished:
   ```
   nextflow run Annotate.nf -w ~/scratch/work_directory -resume
   ```
   
3. If `nextflow` crashes with error `libnet.so: failed to map segment from shared object`, then try to increase the amount of memory in your `salloc` or `sbatch` job.
