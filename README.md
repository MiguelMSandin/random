# A set of random scripts

Random scripts to deal with common fasta files, phylogenetic trees, tab-delimited lists and other fancy stuff. Such as split an alignment in defined positions, concatenate different fasta files, extract/remove sequences, reorder an alignment based on a tree, get the entropy on each position from an alignment, get the consensus sequence from an alignment, ...

Download and move the scripts to you prefered folder (e.g.;```/usr/lobal/bin/```) and start using them (you might have to make them executable: ```chmod +x *.py```).  
  
## A brief description of most useful scripts
  
### For fasta files manipulations ([random/fasta](https://github.com/MiguelMSandin/random/tree/main/fasta))  :  
  
**alignmentConsensus.py**: Creates a consensus of an aligned fasta file.  
**alignmentEntropy.py**: Exports a table with the entropy and other values for every position of an aligned fasta file.  
**fastaConcat.py**: From different fasta files, concatenates the sequences from identical sequence names.  
**fastaRevCom.py**: Exports the reverse complement (or only reversed or only complement) fasta file.  
**fastaSplit.py**: Takes an aligned fasta file and creates several fasta files cut at desired positions.  
**fastaStats.py**: Prints to console certain statistics from a fasta file. Number of sequences, GC content, ...  
**multi2linefasta.py**: From a fasta file where the sequences are saved in different sequences, exports a fasta file with each sequence in one line.  
**multi2line**: Same script as *multi2linefasta.py* but written in bash.  
**sequenceSelect.py**: From a fasta file, selects sequences from a given list or pattern and removes them or extracts them.  
  
### For phylogenetic tree manipulations ([random/phylogenetics](https://github.com/MiguelMSandin/random/tree/main/phylogenetics)):  
  
**fastaReorder.py**: Orders the sequences of a fasta file by border of appearance in an phylogenetic tree.  
**findSeqs.py**: From a phylogenetic tree and a fasta file, finds sequences that do not appear in either file.  
**treeCheckIntruders.py**: From a phylogenetic tree and an attribute file, find sequences resolved within other attribute. Useful to identify badly placed sequences, long branch attraction artifacts and so on.  
**treeColourBranches.py**: Based on a table, colours the branches of a phylogenetic tree, and goes inwards if the colours are monophyletic.  
**treeCountTips.py**: Simply counts the number of tips of one or several phylogenetic trees.  
**treePruneList.py**: Given a list of tips, prunes a phylogenetic tree of those tips.  
**treePruneOutliers.py**: From a phylogenetic tree, finds and prunes branches which length is an outlier relative to the tree.  
**treeRemoveBranchLengths.py**: Simply removes the branch lengths of a phylogenetic tree.  
**treeRoot.py**: Given a list of tips, roots the tree in the last common ancestor of such tips.  
**treeRootOutgroup.py**: Given a list of tips, creates an outgroup in the last common ancestor of such tips, rooting the tree with the outgroup.  
**treeStats.py**: Prints to console certain statistics of a phylogenetic tree. Number of tips, nodes, branch lengths, ...  
**treeTipExtract.py**: Exports the tip names of a phylogenetic tree to a list.  
**treeTipRename.py**: Given a table, changes the tree tip names to the desire output.  

#### Other fancy tools:  
How to cut, trim and concatenate a fasta file composed of two genes [random/wrappers](https://github.com/MiguelMSandin/random/tree/main/wrappers):  
**trimTwoGenes.sh**  
  
Remove reciprocal and identical hits from a table ([random/others](https://github.com/MiguelMSandin/random/blob/main/others/listRemoveReciprocals.py)):  
**listRemoveReciprocals.py**  
  
For example how to create a constraint tree from a fasta file and a newick tree file containing the name of the groups ([random/phylogenetics](https://github.com/MiguelMSandin/random/tree/main/phylogenetics)):  
**checkConstrainTree.py**: From a fasta file and a newick tree file, checks if every tree tip finds its sequences in the fasta file, and exports to a directory the list of sequence names found in the fasta file for every tree tip.  
**checkConstrainTaxa.sh**: Checks each file from the previous directory, and exports a sumarised output. Be aware this script is highly dependent on the structure of your fasta sequence names...  
**buildConstrainTree.sh**: Takes the original tree and the created directory with the sequence names to create the constraint tree to be used for phylogenetic analysis, with each group composed of a polytomy of all the sequences found in the fasta file.  
  
Further details in the help of each script (```script -h```).

Please feel free to report bugs or suggestions. I'm creating these scripts as I need them, so they might also be targeted to my particular needs.  
  
