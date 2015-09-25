# Name: remove_gaps.R
# Auth: u.niazi@imperial.ac.uk
# Date: 24/09/15
# Desc: script to remove gaps from the directory of sequences

source('HPRU_PHE_Flu_Header.R')
cvFiles = dir('Data_external/Reference_sequences/FLU/refset_quasi_bam/')

print(paste('number of files in directory', length(cvFiles)))

# create a data frame of files
lSplit = strsplit(cvFiles, '\\.')
length(lSplit)
df = do.call(rbind, lSplit)
colnames(df) = c('Type', 'Sub_type', 'Strain', 'none')
head(df)

# contingency table
ftable(Type ~ Sub_type, data=df)

dotchart(t(xtabs(~ Type + Sub_type, data=df)))

# repeat the analysis for each subtype
for (st in seq_along(unique(df[,'Sub_type']))){
  csSubtype = unique(df[,'Sub_type'])[st]
  df.sub = matrix(df[df[,'Sub_type'] == csSubtype,], ncol = 4)
  cvFiles.selected = apply(df.sub, 1, function(x) paste(x, collapse = '.'))
  seq = DNAStringSet()
  while(length(cvFiles.selected) > 1000) {
    seq = append(seq, readDNAStringSet(paste('Data_external/Reference_sequences/FLU/refset_quasi_bam/', 
                                             cvFiles.selected[1:1000], sep='')))
    cvFiles.selected = cvFiles.selected[-(1:1000)]
  }
  seq = append(seq, readDNAStringSet(paste('Data_external/Reference_sequences/FLU/refset_quasi_bam/', cvFiles.selected, sep='')))
  seq.2 = sapply(seq_along(seq), function(x) f_XStringReplaceMotif(seq[[x]]))
  seq.2 = DNAStringSet(seq.2)
  names(seq.2) = names(seq)
  dir.create('Results', showWarnings = F)
  dir.create('Results/Seq_to_align', showWarnings = F)
  writeXStringSet(seq.2, filepath = paste('Results/Seq_to_align/', csSubtype, '.fa', sep=''))
  gc()
} # for