# Introduction

This folder contains the main data and results of computational modelling in the paper.

## Directory *trees*:

* trees_all64.pdf:
The plots of all the 64 image trees.

### Directory *4trees*:
This folder contains the seven datasets reported in the main text of the paper. Four of the datasets have image trees.
For each dataset without image tree, there are only two files related to copy numbers.

* dataset_X.txt:
This file contains original copy numbers for each cell in the dataset.


* dataset_X_rvec.txt:
This file contains summary statistics computed from copy numbers.


* dataset_X.nwk:
This file contains original image (phylogenetic) tree, where internal nodes represent cell divisions and terminal nodes represent leaves
The branch length of each branch is in frames (multiply by 3 to get minutes).

* dataset_X_divtree.txt:
The cell lineage tree in the format of edge list, where each node is a cell.
This file is used to derive half tree length.

* dataset_X_divtree_blenratio.txt:
The cell lineage tree in the format of edge list, where each node is a cell.
This file is used to derive the average ratio of the waiting time (branch length) to division of a cell’s parent to that of the cell’s daughter. \
Description of each column:
  1. type: type of current cell division (N: Normal, LC: lagging chromosome, AB: anaphase bridge, MP: multipolar events).
  2. time: cell cycle time of current cell.
  3. after: cell cycle time of one daughter cell.
  4. before: cell cycle time of the parent cell.
  5. after21: cell cycle time of one granddaughter cell.
  6. after22: cell cycle time of the other granddaughter cell.


* dataset_X_full_elist.txt:
This original image (phylogenetic) tree in the format of edge list, with an additional column suggesting the type of a cell division (E: erroneous, N: normal). This file is used in likelihood-ratio test.



#### Data summary
| Dataset in the paper | File name prefix | Comments |
| ----------- | ----------- | ----------- |
| PDTO-9 #1 | "20181127" | 9 errors, 4 LC, 4 AB, 1 MP, 1 death |
| PDTO-9 #4 | "20190409_1" |  13 errors, 5 LC, 7 AB, 1 MP |
| PDTO-9 #5 | "20190409_2" | 7 errors, 3 LC, 4 AB |
| PDTO-9 #7 | "20200304" | 11 errors, 5 LC, 5 AB, 1 MP, 1 death |
| PDTO-9 #2 | "20190310" |
| PDTO-9 #3 | "20190331" |
| PDTO-9 #6 | "20200204" |


### Directory *60trees*:
This folder contains additional 60 trees used in the likelihood-ratio test of birth rate change after mitotic errors in pure birth tree.
For each tree, there are 3 files.

* dataset_X.nwk:
This file contains original image (phylogenetic) tree, where internal nodes represent cell divisions and terminal nodes represent cells.
The branch length of each branch is in frames (multiply by 4.5 to get minutes).

* dataset_X.txt:
The image (phylogenetic) tree in the format of edge list.

* dataset_X_full.txt:
The image (phylogenetic) tree in the format of edge list and an additional column suggesting whether each node (cell division) is normal (N) or has mitotic errors (E).


## Directory *abc_smc*:
This folder contains the results of running ABC SMC on the real data under neutral model (in sudfolder **neutral**) and selection model (in sudfolder **selection**), respectively.


## Directory *dic*:
The folder contains the results of computing DIC (deviance information criterion) on the real data under neutral model and selection model, respectively.


## Directory *power_analysis*:
The folder contains the results of power analysis (model selection under different sets of parameters) on simulated data.


## File *cytoBand_hg19.txt*:
This file is used for converting copy number alterations into arm-level.
It was downloaded from UCSC by command: \
`wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz`.
