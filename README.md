# vep_pipeline

## 1. Set up
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
singularity -B `pwd`/vep_cache:/opt/vep/.vep vep.sif -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep
```
This step may take more than 1h.

### 1.2. Setting up LoFtee

More detailed instructions on how to set up LoFtee are [here](https://github.com/konradjk/loftee) 

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
