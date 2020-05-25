# HGTboot
Program HGTboot implements an iterative DP algorithm for finding acyclic evolutionary scenarios with well-supported and acyclic horizontal gene transfers.
## Requirements:
Python 2.7

## Usage:
python hgtboot.py -s <file with a species tree> -g <file with a gene tree> -g <file with bootstrap trees> [-d .., -l .., -t.., -i .., -h]
  
###### Examples:
For options description:
python hgtboot.py -h 

Running the program with default options:
python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A  

Running the algorithm for 4 iterations:
python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A -i 4

Setting costs for evolutionary events (duplications and losses):
python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A -d 3 -l 0

  
  
