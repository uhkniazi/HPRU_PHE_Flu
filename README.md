# HPRU_PHE_Flu
Seasonal flu sequences analysis

# SCRIPTS

## HPRU_PHE_Flu_Header.R
Header file used with all scripts

## remove_gaps.R
used to remove --- gaps from the set of sequences in the directory. The sequences are named by subtypes, so it uses the 
subtype argument to produce X number of files based on number of subtypes.
NOTE: add command line arguments later

## cluster_alignment_file.R
used to clusters the MSA alignmnet file in fasta format. Imports the file, calculates distance matrix for the sequences
using hamming/edit distance. Displays a tree and PCA plot based on this distance matrix. Reports number of possible clusters
and produces a csv file with each file and assigned cluster label.  
Will break down the clusters further into subclusters and plot the largest subcluster.  
NOTE: add command line and user input facility later.


  
# FUNCTIONS

## f_XStringReplaceMotif
### Desc:
Takes a XString object, finds the motif (usually a gap - by default) and replaces it with nothing (default)
and returns the XString object without the gap motif. Can be used to remove - from long strings.  
### Args:  
oSeq = XString Object, usually DNA or AA string  
motif = '-' Default, can be changed to another motif  
replace = '', the string to replace the motif with  
### Rets:
returns the XString object without the motif

## f_iGetPCAClusterCount
### Desc:
Takes first two components of the PCA and counts possible clusters. The function does this by 
binning the vector of data into X bins, assigning a class label to each bin, counting how many
observations in each bin, and total number of bins with at least one observations. This is calculated
for both the components of the pca matrix, and the max number of bins with at least one observation, in
first or second dimension is reported.
### Args:
pr.out = principal component object returned by prcomp function
### Rets:
returns possible number of clusters in the data