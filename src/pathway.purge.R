setwd('/Volumes/N1/SIEVE')

# Purge by Pathway

############################## DATASETS ########################################
# Variant GQ filtering
GQ.threshold = 20
tmp.singletons<-readRDS('./DATA/snps.combined.info.M2_singletons.rds')

singletons.filtered<-tmp.singletons[which(!is.na(tmp.singletons$GQ) & tmp.singletons$GQ < GQ.threshold & tmp.singletons$generation == 'M2'),]
singletons.filtered.M2<-unique(paste(singletons.filtered$chromosome,singletons.filtered$position,singletons.filtered$REF,singletons.filtered$ALT,sep=':'))

singletons.filtered<-tmp.singletons[which(!is.na(tmp.singletons$GQ) & tmp.singletons$GQ < GQ.threshold),]
singletons.filtered.M2.M5<-unique(paste(singletons.filtered$chromosome,singletons.filtered$position,singletons.filtered$REF,singletons.filtered$ALT,sep=':'))

rm(singletons.filtered)

singletons.M2 <- read.table('./DATA/singletons.M2.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
singletons.M5 <- read.table('./DATA/singletons.M5.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
multix.M5 <- read.table('./DATA/multix.M5.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

#Filter away singletons with scores below threshold
singletons.M2 <- data.frame(singletons.M2, 'variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2, sep=':'), stringsAsFactors = F)
singletons.M5 <- data.frame(singletons.M5, 'variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2, sep=':'), stringsAsFactors = F)

singletons.filtered.M2 <- unique(singletons.filtered.M2,singletons.M2[which(singletons.M2$is.alt == FALSE),]$variant)
singletons.filtered.M2.M5 <- unique(singletons.filtered.M2.M5,c(singletons.M2[which(singletons.M2$is.alt == FALSE),]$variant,singletons.M5[which(singletons.M5$is.ref == FALSE),]$variant))

multix.M5 <- data.frame(multix.M5, 'variant'=paste(multix.M5$chr,multix.M5$pos,multix.M5$allele1,multix.M5$allele2, sep=':'), stringsAsFactors = F)
singletons.M2<-singletons.M2[which(!singletons.M2$variant %in% singletons.filtered.M2.M5),]
singletons.M5<-singletons.M5[which(!singletons.M5$variant %in% singletons.filtered.M2.M5),]
multix.M5<-multix.M5[which(!multix.M5$variant %in% singletons.filtered.M2.M5),]
singletons.M2<-singletons.M2[,1:(ncol(singletons.M2)-2)]
singletons.M5<-singletons.M5[,1:(ncol(singletons.M5)-2)]
multix.M5<-multix.M5[,1:(ncol(multix.M5)-1)]

singletons.M2<-singletons.M2[which((singletons.M2$allele1 == 'G' & singletons.M2$allele2 == 'A') | (singletons.M2$allele1 == 'C' & singletons.M2$allele2 == 'T')),]
singletons.M5<-singletons.M5[which((singletons.M5$allele1 == 'G' & singletons.M5$allele2 == 'A') | (singletons.M5$allele1 == 'C' & singletons.M5$allele2 == 'T')),]
cohort.M2<-sort(unique(c(singletons.M2[!is.na(singletons.M2$host.hom),]$host.hom,singletons.M2[!is.na(singletons.M2$host.het),]$host.het)))
cohort.M5<-sort(unique(c(singletons.M5[!is.na(singletons.M5$host.hom),]$host.hom,singletons.M5[!is.na(singletons.M5$host.het),]$host.het)))

#Cohort
# Make list of selected M2 plants
filtered <- read.table('./DATA/Related_M2_lines.csv', header=F, sep=',',fill=T,stringsAsFactors = F)
filtered <- filtered$V1
fam <- read.table('./DATA/snps.combined.fam', header=F, sep=' ',fill=T,stringsAsFactors = F)
fam<-fam[which(!startsWith(fam$V1,'M5')),]
controls <- fam[which(startsWith(fam$V1,'C')),]$V1
cases <- fam[which(!startsWith(fam$V1,'C')),]$V1
selection <- read.table('./DATA/selection.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
filtered <- unique(c(filtered,sort(unique(selection[which(selection$selected == 0),]$line))))
controls <- setdiff(controls,filtered)
cases <- setdiff(cases,filtered)
cohort.ids <- sort(unique(c(cases)))

# Plant singleton counts M2
tmp.het<-data.frame(table(singletons.M2[which(!is.na(singletons.M2$host.het)),]$host.het), stringsAsFactors = F)
tmp.hom<-data.frame(table(singletons.M2[which(!is.na(singletons.M2$host.hom)),]$host.hom), stringsAsFactors = F)
names(tmp.het)<-c('id','singletons.het')
names(tmp.hom)<-c('id','singletons.hom')
plant.singletons.M2<-merge(tmp.het,tmp.hom, by='id')

# Plant singleton counts M5
tmp.het<-data.frame(table(singletons.M5[which(!is.na(singletons.M5$host.het)),]$host.het), stringsAsFactors = F)
tmp.hom<-data.frame(table(singletons.M5[which(!is.na(singletons.M5$host.hom)),]$host.hom), stringsAsFactors = F)
names(tmp.het)<-c('id','singletons.het')
names(tmp.hom)<-c('id','singletons.hom')
plant.singletons.M5<-merge(tmp.het,tmp.hom, by='id')
rm(tmp.hom,tmp.het)

# Singleton retention
tmp.M2<-data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'),'host.hom'=singletons.M2$host.hom,'host.het'=singletons.M2$host.het)
tmp.M5<-data.frame('variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2,sep=':'),'host.hom'=singletons.M5$host.hom,'host.het'=singletons.M5$host.het)

## Variant chromosome assignment
singleton.chrom<-unique(rbind(data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'),'chromosome'=singletons.M2$chr),data.frame('variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2,sep=':'),'chromosome'=singletons.M5$chr)))

# Variant chromosomal arm assignment
centromeres<-unique(rbind(data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'),'chromosome'=singletons.M2$chr,'pos'=singletons.M2$pos),data.frame('variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2,sep=':'),'chromosome'=singletons.M5$chr,'pos'=singletons.M5$pos)))
centromeres<-merge(centromeres,read.table('./DATA/centromeres.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F), by='chromosome')
centromeres<-data.frame(centromeres[2], 'chromosome.arm'=as.integer(centromeres$pos > centromeres$cen.pos)+1,stringsAsFactors = F)
singleton.chrom<-merge(singleton.chrom,centromeres,by='variant')
rm(centromeres)

singleton.retention<-merge(tmp.M2,tmp.M5, by='variant', all=TRUE)
names(singleton.retention)[2:5] <- c("host.hom.M2","host.het.M2","host.hom.M5","host.het.M5")

tmp.hom<-singletons.M2[which(!is.na(singletons.M2$host.hom)),c(1,2,3,4,5)]
tmp.het<-singletons.M2[which(!is.na(singletons.M2$host.het)),c(1,2,3,4,6)]
names(tmp.hom)[5]<-'host'
names(tmp.het)[5]<-'host'
tmp.M2<-rbind(data.frame(tmp.hom,'genotype.M2'='HOM',stringsAsFactors = F),data.frame(tmp.het,'genotype.M2'='HET',stringsAsFactors = F))

tmp.hom<-singletons.M5[which(!is.na(singletons.M5$host.hom)),c(1,2,3,4,5)]
tmp.het<-singletons.M5[which(!is.na(singletons.M5$host.het)),c(1,2,3,4,6)]
names(tmp.hom)[5]<-'host'
names(tmp.het)[5]<-'host'
tmp.M5<-rbind(data.frame(tmp.hom,'genotype.M5'='HOM',stringsAsFactors = F),data.frame(tmp.het,'genotype.M5'='HET',stringsAsFactors = F))

tmp.M2<-data.frame('variant'=paste(tmp.M2$chr,tmp.M2$pos,tmp.M2$allele1,tmp.M2$allele2,sep=':'),'host'=tmp.M2$host,'genotype.M2'=tmp.M2$genotype.M2)
tmp.M5<-data.frame('variant'=paste(tmp.M5$chr,tmp.M5$pos,tmp.M5$allele1,tmp.M5$allele2,sep=':'),'host'=tmp.M5$host,'genotype.M5'=tmp.M5$genotype.M5)

singleton.retention.mono<-merge(tmp.M2,tmp.M5, by='variant', all=TRUE)
names(singleton.retention.mono)[c(2,4)]<- c("host.M2","host.M5")

purge<-singleton.retention.mono

#Remove variants that go from singleton in M2 to being in multiple plants in M5
purge<-purge[which(!purge$variant %in% paste(multix.M5$chr,multix.M5$pos,multix.M5$allele1,multix.M5$allele2,sep=':')),]

#remove other kinds of outliers such as those that switch host from M2 to M5 etc.
purge<-purge[which(!(!is.na(purge$host.M2) & !is.na(purge$host.M5) & purge$host.M2 != purge$host.M5)),]
purge<-purge[which(!(is.na(purge$host.M2) & !is.na(purge$host.M5))),]
purge<-purge[which(!((!is.na(purge$host.M2) & !purge$host.M2 %in% cohort.ids) | (!is.na(purge$host.M5) & !purge$host.M5 %in% cohort.ids))),]
purge<-purge[c(1,3,5)]

purge[which(purge$genotype.M2 == 'HET'),]$genotype.M2 = '0/1'
purge[which(purge$genotype.M2 == 'HOM'),]$genotype.M2 = '1/1'
purge[which(purge$genotype.M5 == 'HET'),]$genotype.M5 = '0/1'
purge[which(purge$genotype.M5 == 'HOM'),]$genotype.M5 = '1/1'
purge[which(is.na(purge$genotype.M5)),]$genotype.M5 = '0/0'

#Keep only heterozygous variants at M2 generation as the other genotype is fixed and cannot change
purge<-purge[which(purge$genotype.M2 != '1/1'),]
#Remove heterozygous variants at M5 generation sinze their genotype is yet to be fixed to purged or not purged.
purge<-purge[which(purge$genotype.M5 != '0/1'),]

purge$purged = purge$genotype.M5 == '0/0'
purge$purged = as.factor(purge$purged)
purge<-purge[c(1,4)]

# Plant chromosomal arm heterozygous singleton count at M2
tmp<-data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'), 'host'=singletons.M2$host.het, stringsAsFactors=F)
tmp <- merge(tmp,singleton.chrom,by='variant')
tmp<-tmp[which(!is.na(tmp$host)),]

arm1.singletons<-data.frame(stringsAsFactors = F)
arm2.singletons<-data.frame(stringsAsFactors = F)

for (chromosome in unique(tmp$chromosome))
{
  df <- tmp[which(tmp$chromosome == chromosome & tmp$chromosome.arm == 1),]
  df<-data.frame(table(df$host),'chromosome'=chromosome,stringsAsFactors = F)
  names(df)<-c('id','chrom.arm.singletons','chromosome')
  df <- df[c(1,3,2)]
  arm1.singletons<-rbind(arm1.singletons,df)
  
  df <- tmp[which(tmp$chromosome == chromosome & tmp$chromosome.arm == 2),]
  df<-data.frame(table(df$host),'chromosome'=chromosome,stringsAsFactors = F)
  names(df)<-c('id','chrom.arm.singletons','chromosome')
  df <- df[c(1,3,2)]
  arm2.singletons<-rbind(arm2.singletons,df)
}

arm1.singletons<-arm1.singletons[order(arm1.singletons$id, arm1.singletons$chromosome),]
arm2.singletons<-arm2.singletons[order(arm2.singletons$id, arm2.singletons$chromosome),]

arm.singletons<-data.frame(stringsAsFactors = F)
for (chromosome in unique(tmp$chromosome))
{
  arm.singletons <- rbind(arm.singletons,merge(tmp[which(tmp$chromosome.arm == 1 & tmp$chromosome== chromosome),],arm1.singletons[which(arm1.singletons$chromosome == chromosome),],by.x='host',by.y='id')[c(2,5,6)])
  arm.singletons <- rbind(arm.singletons,merge(tmp[which(tmp$chromosome.arm == 2 & tmp$chromosome== chromosome),],arm2.singletons[which(arm2.singletons$chromosome == chromosome),],by.x='host',by.y='id')[c(2,5,6)])
}
names(arm.singletons)[2]<-'chromosome'
rm(tmp,df,arm1.singletons,arm2.singletons)

#PATHWAY
############################## DATASETS ########################################

# Get gene ids for each variant  
annotations <- read.table('./DATA/snps.combined.annotation.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
variants.gene <- data.frame('variant'=paste(annotations$chr,annotations$position,annotations$ref,annotations$alt,sep=':'), 'gene'=annotations$gene, stringsAsFactors=F)

# Remove variants that are not within genes
variants.gene<-variants.gene[which(variants.gene$gene!=""),]

gene.id.translation <- read.table('./DATA/gene.id.translation.tsv', header=F, sep='\t',fill=T,stringsAsFactors = F)
names(gene.id.translation)<-c('BD21.3','BD21')
gene.id.translation$BD21.3<-sapply(strsplit(gene.id.translation$BD21.3,split='.',fixed=T),'[',2)
gene.id.translation$BD21<-toupper(gene.id.translation$BD21)

pathways <- read.table('./DATA/brachypodiumcyc_pathways.20230103', header=T, sep='\t',quote='',fill=T,stringsAsFactors = F)
pathways <- unique(pathways[which(pathways$Gene.id != 'unknown'),c(1,2,7)])

pathways<-merge(pathways,gene.id.translation, by.x='Gene.id', by.y='BD21')[c(2,3,4)]
names(pathways)<-c('pathway.id','pathway.name','gene')
pathways<-unique(pathways)

#----------------------------- TABLE -------------------------------------------
tables<-list()
for (generation in c('M2','M3','M4'))
{
  df <- data.frame()
  for (table in c('1a','1b','2a','2b','3a','3b','4a','4b'))
  {
    tmp <- read.table(paste('./DATA/',generation,'.table.',table,'.tsv',sep=''), header=F, sep='\t',fill=T,stringsAsFactors = F)
    for (y in 1:nrow(tmp))
      for (x in 1:ncol(tmp))
        df<-rbind(df,data.frame('ID'=tmp[y,x],'table'=table,'x'=x,'y'=y,stringsAsFactors = F))
  }
  table <- df
  table<-list(table)
  names(table)<-paste('table',generation,sep='.')
  tables<-c(tables,table)
}

tables<-tables[['table.M3']][c(1,2)]
tables<-tables[which(tables$ID != ""),]

variant.tables<-unique(merge(tables, singleton.retention.mono[c(1,2)], by.x='ID', by.y='host.M2')[c(2,3)])

#----------------------------- COMBINED TABLE ----------------------------------

pathway.purge <- merge(merge(purge, variants.gene, by='variant'),pathways,by='gene',all=T)
pathway.purge[which(is.na(pathway.purge$pathway.id)),]$pathway.id = 'UNKNOWN'
pathway.purge[which(is.na(pathway.purge$pathway.name)),]$pathway.name = 'UNKNOWN'
pathway.purge<-merge(merge(pathway.purge, variant.tables,by='variant'), arm.singletons, by='variant')

############################## ANALYSIS ########################################

results <- list()
res.pathway.purge <- data.frame(stringsAsFactors = F)

for (pathway in unique(pathway.purge$pathway.id))
{
  print(pathway)
  df<-pathway.purge
  df$in.pathway = as.integer(df$pathway.id == pathway)
  genes <- length(unique(df[which(df$pathway.id == pathway),]$gene))
  #res <- glm(formula = in.pathway ~ purged + chromosome+chrom.arm.singletons*table, family = binomial(link = "logit"), data = df)
  res <- glm(formula = purged ~ in.pathway + chromosome+chrom.arm.singletons + table, family = binomial(link = "logit"), data = df)
  results[paste(pathway,'result', sep = '.')] = list(res)
  intercept<-coef(summary(res))[1,1]
  coef<-coef(summary(res))[2,1]
  p<-coef(summary(res))[2,4]
  res.pathway.purge<-rbind(res.pathway.purge, data.frame('pathway'=pathway,'genes'=genes,'intercept'=intercept,'coef'=coef,'p'=p,stringsAsFactors = F))
}

res.pathway.purge<-unique(merge(pathways[c(1,2)],res.pathway.purge, by.x='pathway.id', by.y='pathway'))
res.pathway.purge<-res.pathway.purge[order(res.pathway.purge$p),]
write.table(res.pathway.purge,file='./RESULTS/pathways.purge.results.tsv', sep='\t', quote = F, row.names = F)
saveRDS(results,file='./RESULTS/pathways.purge.results.rds')
View(res.pathway.purge[which(res.pathway.purge$p < 0.05),])

min.pathway.genes <- 11
res.pathway.purge.adjusted<-res.pathway.purge[which(res.pathway.purge$genes >= min.pathway.genes),]
res.pathway.purge.adjusted$p.fdr = p.adjust(res.pathway.purge.adjusted$p, method = 'fdr')
write.table(res.pathway.purge.adjusted,file='./RESULTS/pathways.purge.results.adjusted.tsv', sep='\t', quote = F, row.names = F)
View(res.pathway.purge.adjusted[which(res.pathway.purge.adjusted$p.fdr < 0.05),])

#-----------------------------------------END-----------------------------------
