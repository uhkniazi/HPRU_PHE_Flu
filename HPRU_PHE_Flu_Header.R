# Name: HPRU_PHE_Flu_Header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 24/09/15
# Desc: Header file for libraries and functions



########## libraries
library(Biostrings)




######### global variables
p.old = par()



########## functions
# Function: f_XStringReplaceMotif
# Desc: Takes a XString object, finds the motif (usually a gap - by default) and replaces it with nothing (default)
#       and returns the XString object without the gap motif. Can be used to remove - from long strings
# Args: oSeq = XString Object, usually DNA or AA string
#       motif = '-' Default, can be changed to another motif
#       replace = '', the string to replace the motif with
# Rets: returns the XString object without the motif
f_XStringReplaceMotif = function(oSeq, motif='-', replace=''){
  oMask = maskMotif(oSeq, motif)
  # create views on the masked string
  v = as(oMask, 'Views')
  # invert the mask and use the ranges to get positions for the masked motif
  r = ranges(gaps(v))
  # remove these motifs from the sequence
  return(replaceAt(oSeq, r, replace))
}


# Function: f_iGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns possible number of clusters in the data
f_iGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  return(iClust.count)
}