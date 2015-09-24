# Name: HPRU_PHE_Flu_Header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 24/09/15
# Desc: Header file for libraries and functions



########## libraries
library(Biostrings)





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