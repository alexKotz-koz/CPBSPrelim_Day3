# CPBS Prelims Day 3
> Goal: 

**The assignment for this homework:** _Test_

The outline of this project is as follows:
- Read in a set of next generation sequence (NGS) reads and a query sequence.
- Break each of the NGS sequence reads into k-mers (all possible substrings of length k that are contained in a string).
- Construct a [De Bruijn Graph](https://en.wikipedia.org/wiki/De_Bruijn_graph) using the prefix and suffix of each k-mer. [Concept Overview](https://www.youtube.com/watch?v=TNYZZKrjCSk&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA&index=51). 
- Construct contigs (contiguous sequences) by following all possible paths through the De Brujin Graph.
- Search each contig for the query string.
- Return the longest contig that contains the query string.

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
- Write small description on how to use the program

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
- `--graph`: Show graph or not (optional)

## Requirements
- Python 3.9 or higher. 
    - This project has been tested with python 3.9 thru 3.12.
- Miniconda (see Installation section for further instructions).
- macOS or Linux based operating system.
