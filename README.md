# A set of random scripts

Random scripts to deal with common fasta files, phylogenetic trees, tab-delimited lists and other fancy stuff. Such as split an alignment in defined positions, concatenate different fasta files, extract/remove sequences, reorder an alignment based on a tree, get the entropy on each position from an alignment, get the consensus sequence from an alignment, ...

Download and move the scripts to you prefered folder (e.g.;```/usr/lobal/bin/```) and start using them (you might have to make them executable: ```chmod +x *.py```).  
  
## A brief description of most useful scripts
  
### For fasta files manipulations    
  
**alignmentConsensus.py**: Creates a consensus of an aligned fasta file.  
**alignmentEntropy.py**: Exports a table with the entropy and other values for every position of an aligned fasta file.  
**fastaConcat.py**: From different fasta files, concatenates the sequences from identical sequence names.  
**fastaRevCom.py**: Exports the reverse complement (or only reversed or only complement) fasta file.  
**fastaSplit.py**: Takes an aligned fasta file and creates several fasta files cut at desired positions.  
**fastaStats.py**: Prints to console certain statistics from a fasta file. Number of sequences, GC content, ...  
**multi2linefasta.py**: From a fasta file where the sequences are saved in different sequences, exports a fasta file with each sequence in one line.  
**multi2line**: Same script as *multi2linefasta.py* but written in bash.  
**sequenceSelect.py**: From a fasta file, selects sequences from a given list or pattern and removes them or extracts them.  
  
### For phylogenetic tree manipulations  
  
**fastaReorder.py**: Orders the sequences of a fasta file by border of appearance in an phylogenetic tree.  
**findSeqs.py**: From a phylogenetic tree and a fasta file, finds sequences that do not appear in either file.  
**treeCheckIntruders.py**: From a phylogenetic tree and an attribute file, find sequences resolved within other attribute. Useful to identify badly placed sequences, long branch attraction artifacts and so on.  
**treeColourBranches.py**:  
**treeCountTips.py**:  
**treePruneList.py**:  
**treePruneOutliers.py**:  
**treeRemoveBranchLengths.py**:  
**treeRoot.py**:  
**treeRootOutgroup.py**:  
**treeStats.py**:  
**treeTipExtract.py**:  
**treeTipRename.py**:  

#### More fancy tools:  
For example how to create a constraint tree from a fasta file and a summarized newick tree file:  
**checkConstrainTree.py**:
**checkConstrainTaxa.sh**:
**buildConstrainTree.sh**: 
  
Further details in the help of each script (```script -h```).

Please feel free to report bugs or suggestions.
