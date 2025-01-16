setwd('/Volumes/N1/SIEVE')

# Variant GQ filtering
GQ.threshold = 20
tmp.singletons<-readRDS('./DATA/snps.combined.info.M2_singletons.rds')

singletons.filtered<-tmp.singletons[which(!is.na(tmp.singletons$GQ) & tmp.singletons$GQ < GQ.threshold & tmp.singletons$generation == 'M2'),]
singletons.filtered.M2<-unique(paste(singletons.filtered$chromosome,singletons.filtered$position,singletons.filtered$REF,singletons.filtered$ALT,sep=':'))

singletons.filtered<-tmp.singletons[which(!is.na(tmp.singletons$GQ) & tmp.singletons$GQ < GQ.threshold),]
singletons.filtered.M2.M5<-unique(paste(singletons.filtered$chromosome,singletons.filtered$position,singletons.filtered$REF,singletons.filtered$ALT,sep=':'))

rm(singletons.filtered)

# Singletons
singletons.M2 <- read.table('./DATA//singletons.M2.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
singletons.M5 <- read.table('./DATA/singletons.M5.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
multix.M5 <- read.table('./DATA/multix.M5.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

#Filter away singletons with scores below threshold
singletons.M2 <- data.frame(singletons.M2, 'variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2, sep=':'), stringsAsFactors = F)
singletons.M5 <- data.frame(singletons.M5, 'variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2, sep=':'), stringsAsFactors = F)
multix.M5 <- data.frame(multix.M5, 'variant'=paste(multix.M5$chr,multix.M5$pos,multix.M5$allele1,multix.M5$allele2, sep=':'), stringsAsFactors = F)

singletons.filtered.M2 <- unique(singletons.filtered.M2,singletons.M2[which(singletons.M2$is.alt == FALSE),]$variant)
singletons.filtered.M2.M5 <- unique(singletons.filtered.M2.M5,c(singletons.M2[which(singletons.M2$is.alt == FALSE),]$variant,singletons.M5[which(singletons.M5$is.ref == FALSE),]$variant))

singletons.M2<-singletons.M2[which(!singletons.M2$variant %in% singletons.filtered.M2.M5),]
singletons.M5<-singletons.M5[which(!singletons.M5$variant %in% singletons.filtered.M2.M5),]
multix.M5<-multix.M5[which(!multix.M5$variant %in% singletons.filtered.M2.M5),]
singletons.M2<-singletons.M2[,1:(ncol(singletons.M2)-2)]
singletons.M5<-singletons.M5[,1:(ncol(singletons.M5)-2)]
multix.M5<-multix.M5[,1:(ncol(multix.M5)-1)]

#Keep only GC-tranistion singletons for analysis
singletons.M2<-singletons.M2[which((singletons.M2$allele1 == 'G' & singletons.M2$allele2 == 'A') | (singletons.M2$allele1 == 'C' & singletons.M2$allele2 == 'T')),]
singletons.M5<-singletons.M5[which((singletons.M5$allele1 == 'G' & singletons.M5$allele2 == 'A') | (singletons.M5$allele1 == 'C' & singletons.M5$allele2 == 'T')),]

cohort.M2<-sort(unique(c(singletons.M2[!is.na(singletons.M2$host.hom),]$host.hom,singletons.M2[!is.na(singletons.M2$host.het),]$host.het)))
cohort.M5<-sort(unique(c(singletons.M5[!is.na(singletons.M5$host.hom),]$host.hom,singletons.M5[!is.na(singletons.M5$host.het),]$host.het)))

selection <- read.table('./DATA/selection.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

#Cohort
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
cohort.ids <- sort(unique(cases))

# Keep only GC->AT singletons that are unique to list of selected M2 plants
singletons.hom<-singletons.M2[which(singletons.M2$host.hom %in% cohort.ids),]
singletons.het<-singletons.M2[which(singletons.M2$host.het %in% cohort.ids),]
rm(cases,controls,fam,filtered, selection)

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
tmp.M2<-data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'),'GC.transition'=F,singletons.M2, stringsAsFactors = F)
tmp.M5<-data.frame('variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2,sep=':'),'GC.transition'=F,singletons.M5, stringsAsFactors = F)
tmp.M2[which((tmp.M2$allele1 == 'G' & tmp.M2$allele2 == 'A') | (tmp.M2$allele1 == 'C' & tmp.M2$allele2 == 'T')),]$GC.transition=T
tmp.M5[which((tmp.M5$allele1 == 'G' & tmp.M5$allele2 == 'A') | (tmp.M5$allele1 == 'C' & tmp.M5$allele2 == 'T')),]$GC.transition=T
singleton.transition<-unique(rbind(tmp.M2[c(1,2)], tmp.M5[c(1,2)]))

# Singleton retention
tmp.M2<-data.frame('variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2,sep=':'),'host.hom'=singletons.M2$host.hom,'host.het'=singletons.M2$host.het)
tmp.M5<-data.frame('variant'=paste(singletons.M5$chr,singletons.M5$pos,singletons.M5$allele1,singletons.M5$allele2,sep=':'),'host.hom'=singletons.M5$host.hom,'host.het'=singletons.M5$host.het)

# Variant chromosome assignment
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
purge<-purge[which(purge$genotype.M2 == '0/1'),]

#Remove heterozygous variants at M5 generation since their genotype is yet to be fixed to purged or not purged.
purge<-purge[which(purge$genotype.M5 != '0/1'),]

# Variant annotations
annotations <- read.table('./DATA/snps.combined.annotation.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
annotations <- data.frame('variant'=paste(annotations$chr,annotations$position,annotations$ref,annotations$alt,sep=':'), 'genic'=annotations$genic, 'intergenic'=annotations$intergenic,'cds'=annotations$cds,'utr'=annotations$utr, 'intronic'=F, stringsAsFactors=F)
annotations[which(annotations$genic & !(annotations$cds| annotations$utr)),]$intronic = T

annotations<-rbind(data.frame('variant'=annotations[which(annotations$intergenic == T),]$variant, 'type'='intergenic',stringsAsFactors = F),
                   data.frame('variant'=annotations[which(annotations$intronic == T),]$variant, 'type'='intronic',stringsAsFactors = F),
                   data.frame('variant'=annotations[which(annotations$cds == T),]$variant, 'type'='cds',stringsAsFactors = F),
                   data.frame('variant'=annotations[which(annotations$utr == T),]$variant, 'type'='utr',stringsAsFactors = F))

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
  arm.singletons <- rbind(arm.singletons,merge(tmp[which(tmp$chromosome.arm == 1 & tmp$chromosome== chromosome),],arm1.singletons[which(arm1.singletons$chromosome == chromosome),],by.x='host',by.y='id')[c(2,6)])
  arm.singletons <- rbind(arm.singletons,merge(tmp[which(tmp$chromosome.arm == 2 & tmp$chromosome== chromosome),],arm2.singletons[which(arm2.singletons$chromosome == chromosome),],by.x='host',by.y='id')[c(2,6)])
}
rm(tmp,df,arm1.singletons,arm2.singletons)

#-------------------------- PDS gene distance ----------------------------------
pds.distance <- read.table('./DATA/variants.tss.distance.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
pds.distance<-data.frame('variant'=paste(pds.distance$chromosome,pds.distance$position,pds.distance$allele1,pds.distance$allele2,sep=':'),'tss.distance'=pds.distance$distance,stringsAsFactors = F)
pds.distance$proximal = (pds.distance$tss.distance <= 2000)

############################## VEP SCORES ######################################

#---------------------------- ESM ----------------------------------------------

esm.scores <- read.table('./DATA/esm.scores.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
esm.scores<-data.frame('variant'=paste(esm.scores$chr,esm.scores$pos,esm.scores$allele1,esm.scores$allele2,sep=':'),esm.scores[c(5,6,7)])
esm.scores$score <-round(esm.scores$score,1)
bd.primary.transcript<-read.table('./DATA/BdistachyonBd21_3_537_v1.2.protein_primaryTranscriptOnly.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
esm.scores<-esm.scores[which(esm.scores$transcript %in% bd.primary.transcript$transcript), ]
tally<-data.frame(table('variant'=esm.scores$variant), stringsAsFactors = F)
if (nrow(tally)>0)
  esm.scores<-esm.scores[which(!(esm.scores$variant %in% tally[which(tally$Freq>1), ]$variant)), ]

#----------------------------- PDS ---------------------------------------------

pds.scores <- read.table('./DATA/pds.scores.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
pds.scores$score = round(pds.scores$score,2)
pds.logits <- read.table('./DATA/pds.logits.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

#------------------------------------- SIFT ------------------------------------

sift.scores <- read.table('./DATA/sift.scores.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
bd.primary.transcript<-read.table('./DATA/BdistachyonBd21_3_537_v1.2.protein_primaryTranscriptOnly.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
sift.scores<-sift.scores[which(sift.scores$transcript %in% paste(bd.primary.transcript$transcript, 'v1.2', sep='.')), ]
tally<-data.frame(table('variant'=sift.scores$variant), stringsAsFactors = F)
if (nrow(tally)>0)
  sift.scores<-sift.scores[which(!(sift.scores$variant %in% tally[which(tally$Freq>1), ]$variant)), ]

#---------------------------------- PHYTOEXPR ----------------------------------

phytoxpr.scores <- read.table('./DATA/phytoexpr.scores.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

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
############################## ANALYSIS ########################################

results <- list()
res.purge <- data.frame(stringsAsFactors = F)

#------------------------------- ESM -------------------------------------------
df<-merge(purge, esm.scores, by='variant')
df<-merge(df,singleton.chrom,by='variant')
df<-merge(df,arm.singletons,by='variant')
df<-merge(df,variant.tables,by='variant')
df<-merge(df,singleton.transition, by='variant')
df$purged = df$genotype.M5 == '0/0'
df$purged = as.factor(df$purged)
df$table = as.factor(df$table)
df$chromosome<-paste(df$chromosome,df$chromosome.arm,sep = '.')
df$chromosome = as.factor(df$chromosome)
df$GC.transition = as.factor(df$GC.transition)
res<-glm(df$purged ~ score+chromosome+chrom.arm.singletons+table, data = df, family = "binomial")

results['esm.result'] = list(res)
intercept<-coef(summary(res))[1,1]
coef<-coef(summary(res))[2,1]
p<-coef(summary(res))[2,4]
res.purge<-rbind(res.purge, data.frame('score'='esm',
                                       'intercept'=intercept,
                                       'coef'=coef,
                                       'p'=p,
                                       stringsAsFactors = F))

#------------------------------- SIFT ------------------------------------------
df<-merge(purge, sift.scores, by='variant')
df<-merge(df,singleton.chrom,by='variant')
df<-merge(df,arm.singletons,by='variant')
df<-merge(df,variant.tables,by='variant')
df<-merge(df,singleton.transition, by='variant')
df$purged = df$genotype.M5 == '0/0'
df$purged = as.factor(df$purged)
df$table = as.factor(df$table)
df$chromosome<-paste(df$chromosome,df$chromosome.arm,sep = '.')
df$chromosome = as.factor(df$chromosome)
df$GC.transition = as.factor(df$GC.transition)
res<-glm(df$purged ~ score+chromosome+chrom.arm.singletons*table, data = df, family = "binomial")

results['sift.result'] = list(res)
intercept<-coef(summary(res))[1,1]
coef<-coef(summary(res))[2,1]
p<-coef(summary(res))[2,4]
res.purge<-rbind(res.purge, data.frame('score'='sift',
                                       'intercept'=intercept,
                                       'coef'=coef,
                                       'p'=p,
                                       stringsAsFactors = F))

#------------------------------- PDS -------------------------------------------

df<-merge(purge, pds.logits, by='variant')
df<-merge(df,arm.singletons,by='variant')
df<-merge(df,singleton.chrom,by='variant')
df<-merge(df,variant.tables,by='variant')
df<-merge(df,singleton.transition, by='variant')
df<-merge(df,pds.distance, by='variant')
df$logit.mean.four_tissues = (df$flag_leaf+df$flower+df$panicle+df$young_leaf)/4
df$purged = df$genotype.M5 == '0/0'
df$purged = as.factor(df$purged)
df$table = as.factor(df$table)
df$chromosome<-paste(df$chromosome,df$chromosome.arm,sep = '.')
df$chromosome = as.factor(df$chromosome)
df$GC.transition = as.factor(df$GC.transition)
df$proximal = as.factor(df$proximal)
res<-glm(purged ~ logit.mean.four_tissues+root+chromosome+chrom.arm.singletons+table, data = df, family = "binomial")

results['pds.result'] = list(res)
intercept<-coef(summary(res))[1,1]
coef<-coef(summary(res))[2,1]
p<-coef(summary(res))[2,4]
res.purge<-rbind(res.purge, data.frame('score'='pds.four_tissues',
                                       'intercept'=intercept,
                                       'coef'=coef,
                                       'p'=p,
                                       stringsAsFactors = F))

intercept<-coef(summary(res))[1,1]
coef<-coef(summary(res))[3,1]
p<-coef(summary(res))[3,4]
res.purge<-rbind(res.purge, data.frame('score'='pds.root',
                                       'intercept'=intercept,
                                       'coef'=coef,
                                       'p'=p,
                                       stringsAsFactors = F))

#--------------------------- PHYTOEXPR -----------------------------------------
df<-merge(purge, phytoxpr.scores, by='variant')
df<-merge(df,arm.singletons,by='variant')
df<-merge(df,singleton.chrom,by='variant')
df<-merge(df,variant.tables,by='variant')
df<-merge(df,singleton.transition, by='variant')
df<-merge(df,pds.distance, by='variant')
df$purged = df$genotype.M5 == '0/0'
df$purged = as.factor(df$purged)
df$table = as.factor(df$table)
df$chromosome<-paste(df$chromosome,df$chromosome.arm,sep = '.')
df$chromosome = as.factor(df$chromosome)
df$GC.transition = as.factor(df$GC.transition)
df$proximal = as.factor(df$proximal)
res<-glm(df$purged ~ score+chromosome+chrom.arm.singletons+table, data = df, family = "binomial")

results['phytoEXPR.result'] = list(res)
intercept<-coef(summary(res))[1,1]
coef<-coef(summary(res))[2,1]
p<-coef(summary(res))[2,4]
res.purge<-rbind(res.purge, data.frame('score'='phytoEXPR',
                                       'intercept'=intercept,
                                       'coef'=coef,
                                       'p'=p,
                                       stringsAsFactors = F))

#-------------------------------- SAVE -----------------------------------------

write.table(res.purge,file='./RESULTS/purge.results.tsv', sep='\t', quote = F, row.names = F)
saveRDS(results,file='./RESULTS/purge.results.rds')

#------------------------------- END ----------------------------------------
