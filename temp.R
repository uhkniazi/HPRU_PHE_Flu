# testing script

# test distance distribution
t = sample(m, 10000)
t = sqrt(t)
r = range(t)
s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
r[1] = floor(r[1])
r[2] = ceiling(r[2])
# which distribution can approximate the frequency of reactome terms
hist(t, prob=T, main='degree distribution of type 2 vertices', breaks=s,
     xlab='log degree', ylab='')
# try negative binomial and poisson distributions
# parameterized on the means
dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
dp = dpois(r[1]:r[2], mean(t))
dnr = dnorm(r[1]:r[2], mean(t), sd(t))
lines(r[1]:r[2], dn, col='black', type='b')
lines(r[1]:r[2], dp, col='red', type='b')
lines(r[1]:r[2], dnr, col='blue', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))





library(Biostrings)

#f = file.choose()

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

csSubtype = 'H2N2'
df.h2n2 = df[df[,2] == 'H2N2',]

cvFiles.selected = apply(df.h2n2, 1, function(x) paste(x, collapse = '.'))
seq = readDNAStringSet(paste('Data_external/Reference_sequences/FLU/refset_quasi_bam/', cvFiles.selected, sep=''))
seq.2 = sapply(seq_along(seq), function(x) f_XStringReplaceMotif(seq[[x]]))
seq.2 = DNAStringSet(seq.2)
names(seq.2) = names(seq)
dir.create('Results', showWarnings = F)
dir.create('Results/Seq_to_align', showWarnings = F)

writeXStringSet(seq.2, filepath = paste('Results/Seq_to_align/', csSubtype, '.fa', sep=''))

seq.2 = lapply(seq, f_XStringReplaceMotif(x))


temp = DNAString('ATCG--ATCG')
matchPattern('--', temp)

m = maskMotif(temp, '-')
maskGaps(temp)


# multiple alignment
msa = readDNAMultipleAlignment(file.choose(), format = 'fasta')
sDist = stringDist(as(msa, 'DNAStringSet'), method='hamming')
hc = hclust(sDist, method = 'single')
plot(hc, hang=-1)
m = as.matrix(sDist)
pr.out = prcomp(m, scale=T)

plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2',
main='PCA comp 1 and 2, not normalized')

# how many clusters in data
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

# cut the tree into sub clusters
cut.pt = cutree(hc, k = iClust.count)
# get the label for the largest cluster
iLargest = which.max(as.vector(table(cut.pt)))
df = as.data.frame(table(cut.pt))
iLargest = df[iLargest,1]
ivLargest = which(cut.pt == iLargest)
ivLargest.not = which(cut.pt != iLargest)

# mask these rows in msa
msa.rm = msa
rowmask(msa.rm) = IRanges(start=ivLargest.not, end=ivLargest.not)

# repeat the clustering process
sDist = stringDist(as(msa.rm, 'DNAStringSet'), method='hamming')
hc = hclust(sDist, method = 'single')
plot(hc, hang=-1)
m = as.matrix(sDist)
pr.out = prcomp(m, scale=T)

plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')



## choose appropriate factor
fSamples = as.factor(x.lumi$Day)
fSamples = as.factor(x.lumi$Sample.group.3)
fSamples = as.factor(x.lumi$Study.Group)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)

