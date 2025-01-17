setwd('/Volumes/N1/SIEVE')

############################## DATASETS #######################################

# Variant GQ filtering
GQ.threshold = 20
tmp.singletons<-readRDS('./DATA/snps.combined.info.M2_singletons.rds')

singletons.filtered<-tmp.singletons[which(!is.na(tmp.singletons$GQ) & tmp.singletons$GQ < GQ.threshold & tmp.singletons$generation == 'M2'),]
singletons.filtered.M2<-unique(paste(singletons.filtered$chromosome,singletons.filtered$position,singletons.filtered$REF,singletons.filtered$ALT,sep=':'))

rm(singletons.filtered)

# Get list of M2 singletons GC - > AT
singletons.M2 <- read.table('./DATA/singletons.M2.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)

#Filter away singletons with scores below threshold
singletons.M2 <- data.frame(singletons.M2, 'variant'=paste(singletons.M2$chr,singletons.M2$pos,singletons.M2$allele1,singletons.M2$allele2, sep=':'), stringsAsFactors = F)
singletons.filtered.M2 <- unique(singletons.filtered.M2,singletons.M2[which(singletons.M2$is.alt == FALSE),]$variant)
singletons.M2<-singletons.M2[which(!singletons.M2$variant %in% singletons.filtered.M2),]
singletons.M2<-singletons.M2[,1:(ncol(singletons.M2)-2)]

singletons.M2<-singletons.M2[which((singletons.M2$allele1 == 'G' & singletons.M2$allele2 == 'A') | (singletons.M2$allele1 == 'C' & singletons.M2$allele2 == 'T')),]
cohort.M2<-sort(unique(c(singletons.M2[!is.na(singletons.M2$host.hom),]$host.hom,singletons.M2[!is.na(singletons.M2$host.het),]$host.het)))

selected = TRUE
if (selected)
{
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
  cohort.ids <- sort(unique(c(cases,controls)))
  
  # Keep only GC->AT singletons that are unique to list of selected M2 plants
  singletons.hom<-singletons.M2[which(singletons.M2$host.hom %in% cohort.ids),]
  singletons.het<-singletons.M2[which(singletons.M2$host.het %in% cohort.ids),]
  rm(cases,controls,fam,filtered, selection)
} else {
  singletons.hom<-singletons.M2[which(!is.na(singletons.M2$host.hom)),]
  singletons.het<-singletons.M2[which(!is.na(singletons.M2$host.het)),]
}

singletons.hom<-paste(singletons.hom$chr,singletons.hom$pos,singletons.hom$allele1,singletons.hom$allele2,sep=':')
singletons.het<-paste(singletons.het$chr,singletons.het$pos,singletons.het$allele1,singletons.het$allele2,sep=':')

length(unique(c(singletons.het,singletons.hom)))

# Get gene ids for each variant  
annotations <- read.table('./DATA/snps.combined.annotation.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
variants.gene <- data.frame('variant'=paste(annotations$chr,annotations$position,annotations$ref,annotations$alt,sep=':'), 'gene'=annotations$gene, stringsAsFactors=F)

nrow(variants.gene)
# Remove variants that are not within genes
variants.gene<-variants.gene[which(variants.gene$gene!=""),]

# Get SIFT annotations and keep only BD21.3 transcripts
sift<-read.table('./DATA/BD_SIFTannotations.xls', header=T, sep='\t',fill=T,stringsAsFactors = F)
nrow(sift)
bd.primary.transcript<-read.table('./DATA/BdistachyonBd21_3_537_v1.2.protein_primaryTranscriptOnly.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F)
sift<-sift[which(sift$TRANSCRIPT_ID %in% paste(bd.primary.transcript$transcript, 'v1.2', sep='.')), ]
sift<-data.frame(sift,'variant'=paste(sift$CHROM,sift$POS,sift$REF_ALLELE,sift$ALT_ALLELE,sep=':'),stringsAsFactors = F)
nrow(sift)
# Remove any variant with more than one transcript
tally<-data.frame(table('variant'=sift$variant), stringsAsFactors = F)
sift<-sift[which(!(sift$variant %in% tally[which(tally$Freq>1), ]$variant)), ]
sift<-sift[which((sift$REF_ALLELE == 'G' & sift$ALT_ALLELE == 'A') | (sift$REF_ALLELE == 'C' & sift$ALT_ALLELE == 'T')),]
sift<-data.frame('variant'=sift$variant,'gene'=sift$GENE_ID,'type'=sift$VARIANT_TYPE,stringsAsFactors = F)
nrow(sift)
# Get number of HOM and HET GC->AT missense singleton mutations observed within in each gene
# Merging sift and variants.gene allow us to keep only the relevant singletons, those from M2 and belonging to selected plants etc.
sift<-merge(sift,variants.gene, by='variant')
nrow(sift)
names(sift)[2] <- "gene.id"
names(sift)[4] <- "gene.name"

tmp.het<-data.frame(table('gene'=sift[which(sift$variant %in% singletons.het & sift$type == 'NONSYNONYMOUS'),]$gene.name),stringsAsFactors = F)
names(tmp.het)[2]<-'missense.singletons.het'
tmp.hom<-data.frame(table('gene'=sift[which(sift$variant %in% singletons.hom & sift$type == 'NONSYNONYMOUS'),]$gene.name),stringsAsFactors = F)
names(tmp.hom)[2]<-'missense.singletons.hom'

gene.singletons<-merge(tmp.het,tmp.hom,by='gene',all=T)
gene.singletons[which(is.na(gene.singletons$missense.singletons.het)),]$missense.singletons.het = 0
gene.singletons[which(is.na(gene.singletons$missense.singletons.hom)),]$missense.singletons.hom = 0
gene.singletons$gene=as.character(gene.singletons$gene)

# Tally of all possible GC->AT missense mutations within the CDS region of each gene
gene.GC.missense <- rbind(read.table('./DATA/Bd1.filtered.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F),read.table('./DATA/Bd2.filtered.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F),read.table('./DATA/Bd3.filtered.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F),read.table('./DATA/Bd4.filtered.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F),read.table('./DATA/Bd5.filtered.tsv', header=T, sep='\t',fill=T,stringsAsFactors = F))
gene.GC.missense.tally <- data.frame(table('gene'=gene.GC.missense$gene))
names(gene.GC.missense.tally) <- c('gene','gc.missense')

gene.id.translation <- read.table('./DATA/gene.id.translation.tsv', header=F, sep='\t',fill=T,stringsAsFactors = F)
names(gene.id.translation)<-c('BD21.3','BD21')
gene.id.translation$BD21.3<-sapply(strsplit(gene.id.translation$BD21.3,split='.',fixed=T),'[',2)
gene.id.translation$BD21<-toupper(gene.id.translation$BD21)

pathways <- read.table('./DATA/brachypodiumcyc_pathways.20230103', header=T, sep='\t',quote='',fill=T,stringsAsFactors = F)
pathways <- unique(pathways[which(pathways$Gene.id != 'unknown'),c(1,2,7)])

pathways<-merge(pathways,gene.id.translation, by.x='Gene.id', by.y='BD21')[c(2,3,4)]
names(pathways)<-c('pathway.id','pathway.name','gene')
pathways<-unique(pathways)

gene.missense<-merge(merge(gene.GC.missense.tally,gene.singletons,by='gene',all=T),pathways[c(1,3)],by='gene', all=T)
gene.missense[which(is.na(gene.missense$pathway.id)),]$pathway.id = 'UNKNOWN'
gene.missense[which(is.na(gene.missense$missense.singletons.het)),]$missense.singletons.het = 0
gene.missense[which(is.na(gene.missense$missense.singletons.hom)),]$missense.singletons.hom = 0

############################## ANALYSIS ########################################

#----- ANOVA

results <- list()
res.pathways <- data.frame(stringsAsFactors = F)

for (pathway in unique(gene.missense$pathway.id))
{
  df <- gene.missense
  df$in.pathway = as.integer(df$pathway == pathway)
  genes <- length(unique(df[which(df$pathway == pathway),]$gene))
  df<-data.frame('in.pathway'=df$in.pathway, 'mut.rate'=(df$missense.singletons.hom+df$missense.singletons.het) /df$gc.missense,'mut.rate.het'=df$missense.singletons.het/df$gc.missense,'mut.rate.hom'=df$missense.singletons.hom/df$gc.missense, stringsAsFactors = F)
  
  df$in.pathway = as.factor(df$in.pathway)

  res<-glm(formula = in.pathway ~ mut.rate.hom+mut.rate.het , family = binomial(link = "logit"), data = df)
  res.rand<-glm(formula = in.pathway ~ 1 , family = binomial(link = "logit"), data = df)
  
  res.anova<-anova(res.rand, res, test="Chisq")

  coef.hom=coef(summary(res))[2,1]
  coef.het=coef(summary(res))[3,1]
  p <- res.anova[2,5]

  results[paste(pathway,'result', sep = '.')] = list(res)
  results[paste(pathway,'anova', sep = '.')] = list(res.anova)
  
  res.pathways<-rbind(res.pathways, data.frame('pathway'=pathway,
                                               'genes'=genes,
                                               'coef.hom'=coef.hom,
                                               'coef.het'=coef.het,
                                               'p'=p,
                                               stringsAsFactors = F))
}

res.pathways<-merge(unique(pathways[c(1,2)]), res.pathways, by.x = 'pathway.id', by.y='pathway', all=T)
res.pathways<-res.pathways[order(res.pathways$p),]
write.table(res.pathways,file='./RESULTS/pathways.anova.results.tsv', sep='\t', quote = F, row.names = F)
saveRDS(results,file='./RESULTS/pathways.anova.results.rds')

min.pathway.genes <- 11
res.pathways.adjusted<-res.pathways[which(res.pathways$genes >= min.pathway.genes),c(1,2,3,4,5,6)]
res.pathways.adjusted$p.fdr = p.adjust(res.pathways.adjusted$p, method = 'fdr')

write.table(res.pathways.adjusted,file='./RESULTS/pathways.anova.results.adjusted.tsv', sep='\t', quote = F, row.names = F)
View(res.pathways.adjusted[which(res.pathways.adjusted$p.fdr < 0.05),])

#---- Model 1: Pathway = mut_density, Model 2: mut_density = pathway

results <- list()
res.pathways <- data.frame(stringsAsFactors = F)

for (pathway in unique(gene.missense$pathway.id))
{
  df <- gene.missense
  df$in.pathway = as.integer(df$pathway == pathway)
  genes <- length(unique(df[which(df$pathway == pathway),]$gene))
  df<-data.frame('in.pathway'=df$in.pathway, 'mut.rate'=(df$missense.singletons.hom+df$missense.singletons.het) /df$gc.missense,'mut.rate.het'=df$missense.singletons.het/df$gc.missense,'mut.rate.hom'=df$missense.singletons.hom/df$gc.missense, stringsAsFactors = F)
  
  res3<-lm(formula = mut.rate ~ in.pathway, data = df)
  df$in.pathway = as.factor(df$in.pathway)
  
  res1<-glm(formula = in.pathway ~ mut.rate, family = binomial(link = "logit"), data = df)
  res2<-glm(formula = in.pathway ~ mut.rate.hom+mut.rate.het , family = binomial(link = "logit"), data = df)
  
  if (is.list(res1))
  {
    results[paste(pathway,'result1', sep = '.')] = list(res1)
    coef<-coef(summary(res1))[2,1]
    p<-coef(summary(res1))[2,4]
  } else {
    coef<-NA
    p<-NA
  }
  
  if (is.list(res2))
  {
    results[paste(pathway,'result2', sep = '.')] = list(res2)
    coef.hom=coef(summary(res2))[2,1]
    coef.het=coef(summary(res2))[3,1]
    p.hom=coef(summary(res2))[2,4]
    p.het=coef(summary(res2))[3,4]
  } else {
    coef.hom=NA
    coef.het=NA
    p.hom=NA
    p.het=NA
  }
  
  if (is.list(res3))
  {
    results[paste(pathway,'result3', sep = '.')] = list(res3)
    coef.model2<-coef(summary(res3))[2,1]
    p.model2<-coef(summary(res3))[2,4]
  } else {
    coef.model2<-NA
    p.model2<-NA
  }
  
  res.pathways<-rbind(res.pathways, data.frame('pathway'=pathway,
                                               'genes'=genes,
                                               'coef'=coef,
                                               'p'=p,
                                               'coef.hom'=coef.hom,
                                               'coef.het'=coef.het,
                                               'p.hom'=p.hom,
                                               'p.het'=p.het,
                                               'coef.model2'=coef.model2,
                                               'p.model2'=p.model2,
                                               stringsAsFactors = F))
}

res.pathways<-merge(unique(pathways[c(1,2)]), res.pathways, by.x = 'pathway.id', by.y='pathway', all=T)
res.pathways<-res.pathways[order(res.pathways$p),]
write.table(res.pathways,file='./RESULTS/pathways.results.tsv', sep='\t', quote = F, row.names = F)
saveRDS(results,file='./RESULTS/pathways.results.rds')

min.pathway.genes <- 18
res.pathways.adjusted<-res.pathways[which(res.pathways$genes >= min.pathway.genes),c(1,2,3,4,5)]
res.pathways.adjusted$p.fdr = p.adjust(res.pathways.adjusted$p, method = 'fdr')
write.table(res.pathways.adjusted,file='./RESULTS/pathways.results.adjusted.tsv', sep='\t', quote = F, row.names = F)
View(res.pathways.adjusted[which(res.pathways.adjusted$p.fdr < 0.05),])

#-----------------------------------------END-----------------------------------
