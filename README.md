# Cycads

[![CI](https://github.com/QYanwei/Cycads/actions/workflows/run_test.yml/badge.svg?branch=master)](https://github.com/QYanwei/Cycads/actions/workflows/run_test.yml)

## Description
Cycads long reads quality analyser

## Installation

```
git clone https://github.com/QYanwei/Cycads
cd Cycads
conda env create --file Cycads/environment.yml --name cycads_env
conda activate cycads_env
pip install Cycads
cycads --help
```


## Quick start 

The example below generates HTML report from `test/ecoli.fq.gz`:

  ```
  cycads --fastq test/ecoli.fq.gz --output_dir test/fastq_output
  ```

## Usages

* FASTQ quality control 
  ``` 
  cycads --fastq test/ecoli.fq.gz --output_dir test/fastq_output
  ```
* FASTQ quality control and filtering
  ```
  cycads --fastq test/ecoli.fq.gz --filter --output_dir test/fastq_output
  ```
* FASTQ quality control and alignment-based error analysis
  ```
  cycads --fastq test/ecoli.fq.gz --reference test/ecoli.reference.fasta --output_dir test/alignment_output
  ```
  
* Alignment-based error analysis based on a pre-existing BAM file
  ```
  cycads --bam test/test.bam --output_dir test/bam_output
  ```


## Parameters details


## Example output




