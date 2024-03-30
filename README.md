# Cycads

## Description
Cycads long reads quality analyser

## Installation

```
conda create -n cycads_env -c bioconda -c conda-forge -c defaults jinja2～=3.1.3 matplotlib~=3.8.3 numpy~=1.26.4 pandas~=2.2.1 pyfastx~=2.1.0 pysam~=0.22.0 scipy~=1.12.0 seaborn~=0.13.2 minimap2~=2.17 samtools~=1.11
git clone https://github.com/QYanwei/Cycads

conda activate cycads_env
Cycads/main.py --help
```

## Dependencies

* System requirement:
  
  Lunix/Unix, MacOS
  
  python3.8+ 

* python dependencies: 

  pyfastx==2.1.0 
  
  pysam==0.22.0 
  
  numpy==1.26.4 
  
  pandas==2.2.1 
  
  seaborn=0.13.2 
  
  jinja2=3.1.3 

  scipy

* binary tools: 

  minimap2 (version 2.17-r941) 
  
  samtools (version 1.11,using htslib 1.11) 

## Installation
conda create -n cycads_env -c bioconda -c conda-forge -c defaults jinja2~=3.1.3 matplotlib~=3.8.3 numpy~=1.26.4 pandas~=2.2.1 pyfastx~=2.1.0 pysam~=0.22.0 scipy~=1.12.0 seaborn~=0.13.2 minimap2~=2.17 samtools~=1.11 

git clone https://github.com/QYanwei/Cycads 

conda activate cycads_env 


## Example

* Quick start with a demo dataset
  ```
  cd /path/to/workdir/ 
  python Cycads.py -fq test/ecoli.fq.gz -ref ref/Reference_Ecoli.fasta -name test 
  ```
## Usages

* 1. only fastq quality control 
  ``` 
  python main.py -fq test/ecoli.fq.gz -name test
  ```
* 2. fastq quality control and data filter
  ```
  python main.py -fq test/ecoli.fq.gz -filter -name test
  ```
* 3. fastq quality control and error analysis
  ```
  python main.py -fq test/ecoli.fq.gz -ref ref/Reference_Ecoli.fasta -name test
  ```
   
* 4. fastq quality control and data filter and error analysis
  ```

  python main.py -fq test/ecoli.fq.gz -fitler -ref ref/Reference_Ecoli.fasta -name test

  ```
* 5. only bam error analysis
  ```
  python main.py -bam test/bam -name test 
  ```
## Parameters details


## Demo


## Copyright


