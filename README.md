# Cycads

[![CI](https://github.com/QYanwei/Cycads/actions/workflows/run_test.yml/badge.svg?branch=master)](https://github.com/QYanwei/Cycads/actions/workflows/run_test.yml)

## Description
Cycads long reads quality analyser

## Installation

```
git clone https://github.com/QYanwei/Cycads
conda env create --file Cycads/environment.yml --name cycads_env
conda activate cycads_env
Cycads/main.py --help
```


## Quick start 

The example below generates HTML report from `test/ecoli.fq.gz`:

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


