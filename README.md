BRaVO
====

## Description
Bacterial Relative Variation Outliers (BRaVO) detects features with differential taxonomic abundances, based on the relative variation index among samples.

## Prerequisites
Python packages: numpy, scipy, pandas, scikit-bio, statsmodels
(Test enviroment: Python 2.7)

## Installation
Set up for requirement packages  
`pip install numpy scipy pandas scikit-bio statsmodels`  
  
Clone this repository into your local machine  
`git clone https://github.com/omics-tools/BRaVO.git`  

## Getting Started

**Preparing input files**

**Usage**

`python bravo.py -t count_table.txt -g group_label.txt [-a 0.05] [-p bravo_out] [-o <dirpath>]`

**optional arguments:**

| Flag | Description | Format, Default etc. |
|:-----------|:------------|:------------|
| **-t**       | Taxonomic abundance (count) table (**required**) | comma-delimited table (see sample data)|
| **-g**       | Grouping label for sample (**required**) | comma-delimited table (see sample data) |
| **-a**       | Alpha level at which to control false discoveries | default: 0.05 |
| **-p**       | Prefix for an output file | default: Prefix of an input count table |
| **-o**       | Output directory | default: The directory in the input file |
| **-v**       | show program's version number and exit  | |
| **-h**       | show this help message and exit         | |

## Version

0.0.1 (beta)

## Licence

[MIT] https://github.com/omics-tools/mitoimp/blob/master/LICENSE
