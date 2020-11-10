# Introduction
This program simulates the growth of a patient tumour-derived organoid (PDTO) with a stochastic birth-death branching process to learn mutation rates and selection forces.
<!-- uses ABC (approximate Bayesian computation) approach to fit  -->

It generates cell lineage tree and single-cell copy number profiles (CNPs). It can also simulates bulk CNPs (average CNPs of all cells in an organoid) of multiple organoids.

<!-- Assume there is only cell birth, since there is at most 1 cell death. -->
Only relative arm-/chr-level CNAs are simulated.
All CNAs are assumed to be reciprocal.



# Install
## Required library
* [GSL](https://www.gnu.org/software/gsl/)
* [boost](https://www.boost.org)

## Compilation
Download the source code and then run `make` in folder "code" to generate the simulation program 'simpdo'.



# Usage
Please run `./simpdo -h` for all options.

Option "verbose" controls how much to output.

* When verbose=0, a vector of summary statistics is written to the standard output.
The summary statistics represents
average absolute unique chr-level copy number changes, average absolute unique arm-level copy number changes, half lineage tree length, average ratios of branch length starting from the parent of a node to that starting from the daughter of the node, number of altered arms across all cells (count 1 when one arm is altered in any cell), sum of average CN of all cells across all arms, number of de novo CNAs, and number of arms altered in more than one cell, respectively.
The first four values were used for ABC inference of mutation rates and birth rates (selection coefficients).

* When verbose=1, two files will be output, including the copy numbers for each cell in the final population and summary information of the simulation.

* When verbose=2, five additional files will be output, including the copy numbers, mutations and lineage tree for all cells generated in the branching process.



## Example
### Simulate single-cell CNPs of cells from a single organoid

* Assume neutral evolution:
```
./simpdo -o ./ -e 25 --chr_prob 0.1 --arm_prob 0.1 --birth_rate 0.5 --skip 3 --seed $RANDOM --verbose 0
```


* Assume selection:
```
./simpdo -o ./ -e 25 --model 1 --fitness "1 -0.1"  --chr_prob 0.1 --arm_prob 0.1 --birth_rate 0.5 --skip 3 --seed $RANDOM --verbose 0
```
Here, use "--model 1" to indicate selection is imposed.
Option --fitness "1 -0.1" specifies there is one selection coefficient being -0.1, which indicates negative selection.
By default, gradual selection based on parent cell is assumed.
When a daughter cell accumulates new CNAs after division, its birth rate will be decreased compared with parent cell.


### Simulate bulk CNPs of of multiple organoids

The simulation of multiple organoids is based on parameters inferred from real single-cell CNPs of cells from single organoids.
The output of ABC inference is required.
