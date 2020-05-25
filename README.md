# HGTboot
Program HGTboot implements an iterative DP algorithm for finding acyclic evolutionary scenarios with well-supported and acyclic horizontal gene transfers.
## Requirements:
Python 2.7

## Usage:
python hgtboot.py -s \<file with a species tree> -g \<file with a gene tree> -b \<file with bootstrap trees> [-d \<int>, -l \<int>, -t \<int>, -i \<int>, -h]
  
###### Examples:
<br/>
For options description:

python hgtboot.py -h 

<br/>
Running the program with default options:

python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A  

<br/>
Running the algorithm for 4 iterations:

python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A -i 4

<br/>
Setting costs for evolutionary events (duplication cost = 3, loss cost = 0):

python hgtboot.py -g data/gtree_A -s data/stree -b data/bootstrap_trees_A -d 3 -l 0

  
  
