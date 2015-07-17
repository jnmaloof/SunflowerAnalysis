setwd("~/Documents/Lab Notebook support/2015/SunflowerAnalysis")

# load the dat
counts <- read.delim("~/Dropbox/RNA-Seq Data/HAB1P1_All6Lanes_HAB2P1_All7Lanes_raw_data_merged_mod.txt",row.names=1)
dim(counts)
head(counts)

# get rid of the unmatched reads
counts<-counts[rownames(counts)!="*",]

#check the data
head(counts)
counts[is.na(counts)] <- 0 # NA should be 0
names(counts)
summary(counts)
sort(colSums(counts))
hist(colSums(counts))
barplot(sort(colSums(counts)))
counts <- counts[,colSums(counts) > 3000000] # remove samples with less than 3M reads.  In this case the ones with unmatched barcodes
head(counts)
names(counts)
#te <- unlist(strsplit(names(counts),"_"))[2]

#assign sample descriptors
samples <- data.frame(
  file = names(counts),
  side = sub("(HA2015_)(E|W)([[:digit:]]+)([A-C])(F|S+)(.bam)", "\\2", names(counts)),  # Two sides represented as East and West
  time = sub("(HA2015_)(E|W)([[:digit:]]+)([A-C])(F|S+)(.bam)", "\\3", names(counts)),  # 12 time-points represented by numbers 1-12
  rep = sub("(HA2015_)(E|W)([[:digit:]]+)([A-C])(F|S+)(.bam)", "\\4", names(counts)), # Three replicates represented by A, B, and C
  day = sub("(HA2015_)(E|W)([[:digit:]]+)([A-C])(F|S+)(.bam)", "\\5", names(counts)) # 2 days represented by F(first day) and S(second day))
)
samples

samples$group <- paste(samples$side,samples$time,samples$day,sep="_")

#require more than 10 coutns in more than 6 samples
counts2 <- counts[rowSums(counts > 10) >= 6,]

#start the analysis; normalize the data
library(edgeR)
dge.data <- DGEList(counts=counts2, group=samples$group)
dim(dge.data) 
dge.data <- calcNormFactors(dge.data, method = "TMM")
dge.data$samples # look at the normalization factors

dim(dge.data)

#plotMDS(dge.data, method = "bcv")

#our model for the experiment
design <- model.matrix(~side*time+day,data = samples)
rownames(design) <- samples$file
design

#calculate dispersions
#First the overall dispersion
dge.data <- estimateGLMCommonDisp(dge.data,design,verbose = TRUE)

#Then a trended dispersion based on count level
dge.data <- estimateGLMTrendedDisp(dge.data,design)

#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
dge.data <- estimateGLMTagwiseDisp(dge.data,design)

#We can examine this with a plot
plotBCV(dge.data)

save.image(file="sunflower.Rdata")

#Fit the full model to the data
fit <- glmFit(dge.data, design)

save.image(file="sunflower.Rdata")

# a list to hold the results of separate time points
lrt.results <- list()

# Here I will drop each "side" by "time" interaction one by one.
for(time in grep(":",colnames(design),value=TRUE)) {
  lrt.results[[make.names(paste(time,"lrt"))]] <- glmLRT(fit,coef=time)
}

names(lrt.results)

sapply(lrt.results,function(x) {
  summary(decideTestsDGE(x,p.value=0.05))
})

# what about dropping all interaction coefficients at once?  This should be more sensitive than the former analysis.
lrt.all.interactions <- glmLRT(fit,coef=grep(":",colnames(design)))

lrt.all.results <- topTags(lrt.all.interactions,n=Inf)$table[topTags(lrt.all.interactions,n=Inf)$table$FDR<.05,]
# 48 genes

# now lets plot them
genes.of.interest <- log2(cpm(dge.data)[row.names(lrt.all.results),]+1)

head(genes.of.interest)

dim(genes.of.interest)

library(reshape2)
library(ggplot2)

genes.of.interest.melt <- melt(genes.of.interest)

head(genes.of.interest.melt)

genes.of.interest.melt <- merge(genes.of.interest.melt,samples,by.x="Var2",by.y="file")
summary(genes.of.interest.melt)
genes.of.interest.melt$time <- as.numeric(as.character(genes.of.interest.melt$time))

head(genes.of.interest.melt)

pdf("genesOfInterest.pdf",width=10,height=8)
for(gene in unique(genes.of.interest.melt$Var1)) {
  pl <- ggplot(genes.of.interest.melt[genes.of.interest.melt$Var1==gene,],aes(x=time,y=value,color=side,fill=side))
  pl <- pl + geom_smooth()
  pl <- pl + geom_point()
  pl <- pl + ggtitle(gene) + ylab("log2(cpm+1)")
  print(pl)
}
dev.off()


## Repeat but include day:time interaction.  This should help account for day-to-day variance

design2 <- model.matrix(~side*time+day*time,data = samples)
rownames(design2) <- samples$file
design2
colnames(design2)
design2 <- design2[,colnames(design2)!="time10:dayS"] #necessary because we are missing this time/day combo

#First the overall dispersion
dge.data2 <- estimateGLMCommonDisp(dge.data,design2,verbose = TRUE)

#Then a trended dispersion based on count level
dge.data2 <- estimateGLMTrendedDisp(dge.data2,design2)

#And lastly we calculate the gene-wise dispersion, using the prior estimates to "squeeze" the dispersion towards the common dispersion.
dge.data2 <- estimateGLMTagwiseDisp(dge.data2,design2)

#We can examine this with a plot
plotBCV(dge.data2)

fit2 <- glmFit(dge.data2, design2)

save.image(file="sunflower.Rdata")

lrt.results2 <- list()

for(time in grep(":time",colnames(design2),value=TRUE)) {
  print(time)
  lrt.results2[[make.names(paste(time,"lrt"))]] <- glmLRT(fit2,coef=time)
}

names(lrt.results2)

sapply(lrt.results2,function(x) {
  summary(decideTestsDGE(x,p.value=0.05))
})

# what about dropping all interaction coefficients

lrt.all.interactions2 <- glmLRT(fit2,coef=grep(":time",colnames(design2)))

lrt.all.results2 <- topTags(lrt.all.interactions2,n=Inf)$table[topTags(lrt.all.interactions2,n=Inf)$table$FDR<.05,]
# 76 genes

# now lets plot them

genes.of.interest2 <- log2(cpm(dge.data2)[row.names(lrt.all.results2),]+1)

head(genes.of.interest2)

dim(genes.of.interest2)

library(reshape2)
library(ggplot2)

genes.of.interest.melt2 <- melt(genes.of.interest2)

head(genes.of.interest.melt2)

genes.of.interest.melt2 <- merge(genes.of.interest.melt2,samples,by.x="Var2",by.y="file")
summary(genes.of.interest.melt2)
genes.of.interest.melt2$time <- as.numeric(as.character(genes.of.interest.melt2$time))

head(genes.of.interest.melt2)

pdf("genesOfInterest2.pdf",width=10,height=8)
for(gene in unique(genes.of.interest.melt2$Var1)) {
  pl <- ggplot(genes.of.interest.melt2[genes.of.interest.melt2$Var1==gene,],aes(x=time,y=value,color=side,fill=side))
  pl <- pl + geom_smooth()
  pl <- pl + geom_point()
  pl <- pl + ggtitle(gene) + ylab("log2(cpm+1)")
  print(pl)
}
dev.off()


