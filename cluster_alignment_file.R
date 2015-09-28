# Name: cluster_alignment_file.R
# Auth: u.niazi@imperial.ac.uk
# Date: 24/09/15
# Desc: script to cluster an MSA file

source('HPRU_PHE_Flu_Header.R')

# multiple alignment
csFile = file.choose() # replace this with command line input in future
# read the file
msa = readDNAMultipleAlignment(csFile, format = 'fasta')
# cluster the data using hamming/edit distance
sDist = stringDist(as(msa, 'DNAStringSet'), method='hamming')

# make copy of variable to use at end
sDist.bk = sDist

hc = hclust(sDist, method = 'single')
plot(hc, hang=-1, label=F, sub='', xlab='')

# calculate number of clusters using PCA
m = as.matrix(sDist)
pr.out = prcomp(m, scale=T)

plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')

# get number of possible clusters
iClust.count = f_iGetPCAClusterCount(pr.out)

print(paste('Number of possible clusters = ', iClust.count))

# cut the tree into sub clusters
cut.pt = cutree(hc, k = iClust.count)
dfReport = data.frame(Seq_name=names(cut.pt), Cluster=cut.pt)
csFile.report = paste(csFile, '.clusters.csv', sep='')
write.csv(dfReport, file=csFile.report)

print(paste('Cluster Distribution'))
print(as.data.frame(table(cut.pt)))

# get the label for the largest cluster
iLargest = which.max(as.vector(table(cut.pt)))
df = as.data.frame(table(cut.pt))
iLargest = df[iLargest,1]
ivLargest = which(cut.pt == iLargest)
ivLargest.not = which(cut.pt != iLargest)

print(paste('Number of Sequences in Largest Cluster = ', length(ivLargest)))
ivLargest.last = 1
#### repeat clustering but on largest clusters to make subclusters
while(TRUE){
  # mask these rows in msa
  msa.rm = msa
  rowmask(msa.rm) = IRanges(start=ivLargest.not, end=ivLargest.not)
  msa = msa.rm
  
  # repeat the clustering process
  sDist = stringDist(as(msa.rm, 'DNAStringSet'), method='hamming')
  hc = hclust(sDist, method = 'single')
  #plot(hc, hang=-1, label=F, sub='', xlab='')
  m = as.matrix(sDist)
  pr.out = prcomp(m, scale=T)
  
  #plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2', main='PCA comp 1 and 2')
  
  # get number of possible clusters
  iClust.count = f_iGetPCAClusterCount(pr.out)
  print(paste('Number of possible clusters = ', iClust.count))
  
  # cut the tree into sub clusters
  cut.pt = cutree(hc, k = iClust.count)
    
  # get the label for the largest cluster
  iLargest = which.max(as.vector(table(cut.pt)))
  df = as.data.frame(table(cut.pt))
  iLargest = df[iLargest,1]
  ivLargest = which(cut.pt == iLargest)
  ivLargest.not = which(cut.pt != iLargest)
  
  print(paste('Number of Sequences in Largest Cluster = ', length(ivLargest)))
  # check if clusters can not be made any smaller
  if (length(ivLargest) == length(ivLargest.last)) break
  ivLargest.last = ivLargest
}

# get the cutpoint for the tree when it has been split into minimum number of clusters
cut.ht = max(hc$height)

sDist = sDist.bk
hc = hclust(sDist, method = 'single')

cut.pt = cutree(hc, h = cut.ht)
# maximum possible number of clusters
print(paste('Maximum possible clusters = ', length(unique(cut.pt))))
dfReport = data.frame(Seq_name=names(cut.pt), Cluster=cut.pt)
csFile.report = paste(csFile, '.clusters_max.csv', sep='')
write.csv(dfReport, file=csFile.report)

