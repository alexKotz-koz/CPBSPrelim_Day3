# CPBS Prelims Day 3
> Metagenomic Sequence Assembler and Virome Characterizer 

**Goal:** Develop a general method for characterizing the viral communities in a metagenomic sample. The input to your method is a short-read fastq file (examples below) containing sequences from a metagenomic sample. 

The outline of this project is as follows:
1. Read in the metagenomic biosample (.fastq) and known viruses (.fasta) 
2. Perform Quality Control on the metagenomic sample:
    - Convert existing quality control ascii values to a Phred score and perform pruning of low quality reads
3. Break the biosample reads into k-mers (user defined)
4. Create a de Bruijn Graph from the pool of biosample k-mers
5. Create contigs via a Depth-First Search traverse across the de Bruijn Graph
6. Search the contigs for substrings of the viral sequences
    - Three implementations exist: 

        a. Standard search string strategy (searchForViruses.py)

        b. Leverage a Smith-Waterman local alignment algorithm for aligning each contig against each virus (searchForViruses_SW.py) ~Non-functional

        c. Improve the Smith-Waterman algorithm by implementing multiprocessing (via the concurrent.futures module) (searchForViruses_SW_PP.py) ~ Non-functional

## Installation

OS X & Linux:
1. Clone or download the repository.
2. Set up miniconda environment:
    - If miniconda is not installed on the local machine, please follow the steps outlined here before continuing: [Miniconda installation](https://docs.anaconda.com/free/miniconda/)
    - Once miniconda is installed, create the conda environment by copying this command into a shell (terminal) with an active base conda environment:
        ```sh
        conda env create -f conda_env.yml
        ```
    - Then activate the new conda environment:
        ```sh
        conda activate conda_env
        ```

## Usage example
- 

1. Open a terminal and navigate to /src/:
```sh
cd src
```
2. Use the following command to run the project:
```sh
python3 main.py [options]
```
Options:
- `-h, --help`: Show help menu
- `-biosample`: Metagenomic biosample file (required)
- `-k`: User defined size of the k-mers (required)

## Requirements
- Python 3.9 or higher. 
    - This project has been tested with python 3.9 thru 3.12.
- Miniconda (see Installation section for further instructions).
- macOS or Linux based operating system.
