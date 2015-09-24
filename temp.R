# testing script
library(Biostrings)

f = file.choose()

# check individual files in directory
cvFiles = dir('Data_external/Reference_sequences/FLU/refset_quasi_bam/')
length(cvFiles)
head(cvFiles)

lSplit = strsplit(cvFiles, '\\.')
length(lSplit)
df = do.call(rbind, lSplit)
head(df)
ftable(V2 ~ V1, data=df)
ftable(V1 ~ V2, data=df)
dotchart(t(xtabs(~ V1 + V2, data=df)))

by(df, df[,2], summary)

df.h2n2 = df[df[,2] == 'H2N2',]

cvFiles.selected = apply(df.h2n2, 1, function(x) paste(x, collapse = '.'))
seq = readDNAStringSet(paste('Data_external/Reference_sequences/FLU/refset_quasi_bam/', cvFiles.selected, sep=''))


temp = DNAString('ATCG--ATCG')
matchPattern('--', temp)

m = maskMotif(temp, '-')
maskGaps(temp)
