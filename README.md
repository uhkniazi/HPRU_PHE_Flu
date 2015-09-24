# HPRU_PHE_Flu
Seasonal flu sequences analysis


  
  
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

