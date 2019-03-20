BRaVO
====
Bacterial Relative Variation Outliers (BRaVO) detects features with differential taxonomic abundances, based on the relative variation index among samples.


## Requirement packages
numpy, scipy, pandas, scikit-bio, statsmodels (Test enviroment: Python 2.7)

## Installation
Clone this repository into your local machine  
`git clone https://github.com/omics-tools/BRaVO.git`  

Set up for requirement packages  
`pip install numpy scipy pandas scikit-bio statsmodels`  

## Getting started

### Preparing input files

Please refer to [example data](https://github.com/omics-tools/BRaVO/examples
) for details of file format.

#### Taxonomic abundance table (comma-separated)

|  | ***Sample 1*** | ***Sample 2*** | ・・・ | ***Sample j*** |
|:-----------:|:------------:|:------------:|:------------:|:------------:|
| ***Species 1 (or OTU 1)***      |*X<sub>11</sub>* |*X<sub>12</sub>*|・・・|*X<sub>1j</sub>*|
| ***Species 2 (or OTU 2)***      |*X<sub>21</sub>*|*X<sub>22</sub>*|・・・ |*X<sub>2j</sub>*|
| **・** <br> **・**<br>**・**<br>|・<br>・<br>・<br>|・<br>・<br>・<br>|・&emsp;&emsp;<br>・<br>&emsp;&emsp;・<br>|・<br>・<br>・<br>|
| ***Species i (or OTU i)*** |*X<sub>i1</sub>*|*X<sub>i2</sub>*|・・・|*X<sub>ij</sub>*|

Note: *X* is the number of counts for a species.

#### Group label table (comma-separated)

|***Group*** | ***Sample*** | 
|:-----------:|:------------:|
| A | Sample 1 |
| A | Sample 2 |
| ・ <br> ・<br>・<br>|・<br>・<br>・<br>|
| B | Sample j-1 |
| B | Sample j |

Note: Changing the header names of "Group" and "Sample" may cause a format error.

### Usage

`python bravo.py -t count_table.txt -g group_label.txt [-a 0.05] [-p bravo_out] [-o <dirpath>]`

**optional arguments:**

| Flag | Description | Format, Default etc. |
|:-----------|:------------|:------------|
| **-t**       | Taxonomic abundance (count) table (**required**) | comma-delimited csv (see sample data)|
| **-g**       | Grouping label for sample (**required**) | comma-delimited table (see sample data) |
| **-a**       | Alpha level at which to control false discoveries | default: 0.05 |
| **-p**       | Prefix for an output file | default: Prefix of an input count table |
| **-o**       | Output directory path | default: The directory in the input file |
| **-v**       | Program's version number and exit  | |
| **-h**       | Help message and exit         | |

## Version

0.0.1 (beta)

## Licence

[MIT] https://github.com/omics-tools/mitoimp/blob/master/LICENSE
