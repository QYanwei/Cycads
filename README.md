# Cycads

## Description
Cycads long reads quality analyser

## Dependencies

* system requirement: 
  Lunix/Unix, MacOS
  python3.8+ 

* python dependencies: 

  pyfastx==2.1.0 
  
  pysam==0.22.0 
  
  numpy==1.26.4 
  
  pandas==2.2.1 
  
  seaborn=0.13.2 
  
  jinja2=3.1.3 

* third-party tools: 

  minimap2 (version 2.17-r941) 
  
  samtools (version 1.11,using htslib 1.11) 
  
  pyfastx (version 2.1.0) 

(tips: Please alias the full path of those tools into the Cycads/tool/ folder, which can be found by Cycads easily.) 


## Example

* quickly use
  ```
  cd /path/to/workdir/ 
  python Cycads.py -fq test/ecoli.fq.gz -ref ref/Reference_Ecoli.fasta -name test 
  ```
## Usages

* 1. only fastq quality control 
  ``` 
  python Cycads.py -fq test/ecoli.fq.gz -name test
  ```
* 2. fastq quality control and data filter
  ```
  python Cycads.py -fq test/ecoli.fq.gz -filter -name test
  ```
* 3. fastq quality control and error analysis
  ```
  python Cycads.py -fq test/ecoli.fq.gz -ref ref/Reference_Ecoli.fasta -name test
  ```
   
* 4. fastq quality control and data filter and error analysis
  ```
  python Cycads.py -fq test/ecoli.fq.gz -fitler -ref ref/Reference_Ecoli.fasta -name test
  ```
* 5. only bam error analysis
  ```
  python Cycads.py -bam test/bam -name test 
  ```
## Parameters details


## Demo


## Copyright


