# ANALYSIS OF DIFFERENTIALLY REGULATED GENES IN POTATO SAMPLES

library(edgeR)  # load library

rm(list=ls()) # clear workspace

# LOAD DATA from featureCounts output

fc <- read.table('/Users/igorfesenko/Solanum_phytozome_read_counts', sep='\t', header=TRUE)

colnames(fc)
dim(fc)

# change names of columns

colnames(fc) <-  sub('.bam','',colnames(fc))
colnames(fc) <-  sub('bam.','',colnames(fc))

colnames(fc)

#keep only counts data

filt_fc <-  fc[,c("Geneid", "rnaseq10_s","rnaseq11_s","rnaseq12_s","rnaseq13_s","rnaseq14_s",
                  "rnaseq15_s","rnaseq16_s","rnaseq1_s","rnaseq2_s","rnaseq3_s","rnaseq4_s","rnaseq5_s",
                  "rnaseq6_s","rnaseq7_s","rnaseq8_s","rnaseq9_s")]
#change columns order
colorder <-  c("Geneid","rnaseq1_s","rnaseq2_s","rnaseq3_s","rnaseq4_s","rnaseq5_s",
               "rnaseq6_s","rnaseq7_s","rnaseq8_s","rnaseq9_s","rnaseq10_s","rnaseq11_s","rnaseq12_s","rnaseq13_s","rnaseq14_s",
               "rnaseq15_s","rnaseq16_s")

filt_fc_col <- filt_fc[, colorder]

#change row names to Gene ID
row.names(filt_fc_col) <- filt_fc_col$Geneid
filt_fc_col$Geneid <- NULL

head(filt_fc_col)

# save data to file
saveRDS(filt_fc_col,'counts.Rdata')

#Creating matrix with info about samples

meta <-  data.frame(infection=c("mock","mock","mock","mock","mock","mock","mock","mock","PVY","PVY","PVY","PVY","PVY","PVY","PVY","PVY"),
                    temp=c(22,22,22,22,28,28,28,28,22,22,22,22,28,28,28,28),
                    col=c("green","green","green","green","orange","orange","orange","orange","blue","blue","blue","blue","red","red","red","red"))
rownames(meta) = colnames(filt_fc_col)

head(meta)

# FILTRATION BASED ON COVERAGE

# distribution according to reads coverage
hist(log10(1+apply(filt_fc_col,1,mean)),20)
abline(v=log10(1+10),col='red') # threshold in 10 reads

# genes without expression
table(apply(filt_fc_col,1,mean) == 0)

# genes having coverage above 10 reads
table(apply(filt_fc_col,1,mean) >= 10) 

# filtering genes  - above 10 reads per sample
fc10 = filt_fc_col[apply(filt_fc_col,1,mean) >= 10,]
table(apply(fc10,1,mean) >= 10)

# save data
saveRDS(fc10,'./RNAseq/Rdata/filtered_counts.Rdata')
saveRDS(meta,'./RNAseq/Rdata/meta.Rdata')


# DE GENES ANALYSIS

#correlation between samples
cors = cor(fc10,m='sp') # Spearman correlation
corp = cor(fc10,m='p') # Pearson 
corlp = cor(log(1+fc10),m='p') # log Pearson

# Multidimensional scaling

mdss =  cmdscale(1-cors ,k = 2) # Spearmen 
mdsp =  cmdscale(1-corp ,k = 2) # Pearson
mdslp = cmdscale(1-corlp,k = 2) # log Pearson

par(mfrow=c(2,2)) # visualisation of MDS

plot(mdss ,pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2', main='Spearman', bty='n')
plot(mdsp, pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='Pearson'    ,bty='n')
plot(mdslp, pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='log Pearson',bty='n')

# Using edgeR for analysis

counts = DGEList(fc10) # create edgeR object
counts$samples

#filtering of low-expressed genes based CPM-calculated values
keep <- rowSums(cpm(counts) > 0.5) >= 5
table(keep)

filt_counts <- counts[keep, , keep.lib.sizes = FALSE]
dim(filt_counts)

# RLE normalisation
edger = calcNormFactors(filt_counts,method='RLE')
edger$samples

# get and out RLE-norm CPMs
cpm = cpm(edger) # cpm
write.csv(cpm, file = "./RNAseq/cpm_values.csv", row.names = TRUE)
head(cpm)

#Creating matrix design
# we are going to find differences between infected and control plants as well as
#temperature (temp) and complex stress conditions (infection:temp)
design = model.matrix(~ infection + temp + infection:temp,data = meta)
design

# biological coefficient of variation (BCV)
edger = estimateDisp(edger,design)
plotBCV(edger)
edger$common.dispersion

# fit GLM model
glm <-  glmQLFit(edger,design)
plotQLDisp(glm)

# calculating DE genes during infection
pv_inferction <-  glmLRT(glm,2)
#summary statistics
summary(decideTests(pv_inferction))

topTags(pv_inferction, n=20) #most diff genes

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(pv_inferction)
abline(h=c(-1, 1), col="blue")

# calculate FDR
pv_inferction$table$FDR <- p.adjust(pv_inferction$table$PValue, method="BH")
pv_inferction$table[1:5,]
sum(pv_inferction$table$FDR<0.05) #the number of DE genes

# p-value distribution
hist(pv_inferction$table$PValue)

#---------------------------------------
# calculating DE genes during heat shock
pv_temp = glmLRT(glm,3)
#summary statistics
summary(decideTests(pv_temp))

topTags(pv_temp, n=20) #наиболее меняющиеся гены

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(pv_inferction)
abline(h=c(-1, 1), col="blue")

# calculate FDR
pv_temp$table$FDR <- p.adjust(pv_temp$table$PValue, method="BH")
sum(pv_temp$table$FDR<0.05) #the number of DE genes

pv_temp$table[1:5,]

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(pv_temp)
abline(h=c(-1, 1), col="blue")

#p-value distribution
hist(pv_temp$table$PValue)

#------------------------------------------

# DE genes under complex stress conditions
pv_complex = glmLRT(glm,4)
#summary statistics
summary(decideTests(pv_complex))

topTags(pv_complex, n=20)

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(pv_complex)
abline(h=c(-1, 1), col="blue")

# calculating FDR
pv_complex$table$FDR <- p.adjust(pv_complex$table$PValue, method="BH")
sum(pv_temp$table$FDR<0.05) #the number of DE genes

pv_complex$table[1:5,]

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(pv_complex)
abline(h=c(-1, 1), col="blue")

# p-value distibution
hist(pv_complex$table$PValue)

#CPM top DE genes in complex stress conditions
top <- rownames(topTags(pv_complex, n=20))
cpm(edger)[top,]

#OUT TABLES
write.csv(pv_complex$table, file = "/Volumes/GoogleDrive/Мой диск/solanum/DIFF/complex.csv", row.names = TRUE)

