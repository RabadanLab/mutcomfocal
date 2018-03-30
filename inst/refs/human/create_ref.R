#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
ref_name <- args[1]

message(paste("Working with reference", ref_name))

# load libraries
# build reference library as shown in https://bioinformatics.stackexchange.com/a/3951/1771
if (ref_name == "hg19") {
    library(Homo.sapiens.hg19)
    ref_lib <- Homo.sapiens.hg19
} else if (ref_name == "hg38") {
    library(Homo.sapiens.hg38)
    ref_lib <- Homo.sapiens.hg38
}
# version 2.10.1
library(plyr)

# get location of all genes
symbols <- keys(ref_lib, keytype="SYMBOL")
gene <- select(ref_lib, keys=symbols, cols <-c('SYMBOL', 'TXCHROM', 'TXSTART', 'TXEND'), keytype='SYMBOL')
colnames(gene) <- c("SYMBOL", "CHRLOCCHR", "CHRLOC", "CHRLOCEND")
# remove chr prefix
gene$CHRLOCCHR <- substring(gene$CHRLOCCHR, 4)

rm(symbols)

# if there are many locations listed for a gene, take the first one
gene <- gene[!duplicated(gene["SYMBOL"]),]

# remove genes that are not located
gene <- na.omit(gene)

# positivize the positions of reverse strand genes
gene$CHRLOC <- abs(gene$CHRLOC)
gene$CHRLOCEND <- abs(gene$CHRLOCEND)

# keep only genes on the 22 autosomes
gene <- gene[which(gene$CHRLOCCHR %in% 1:22),]

# rename columns
gene <- with(gene, data.frame(gene=SYMBOL, chr=CHRLOCCHR, begin=CHRLOC, end=CHRLOCEND))

# number of genes per chromosome
chr <- ddply(gene, .(chr), nrow)
names(chr) <- c("idx", "num_genes")
chr$name <- chr$idx

# get centromere positions
# by combining p11.1 and q11.1 from cytoBand
centromere_text <- system(paste("curl http://hgdownload.cse.ucsc.edu/goldenPath/", ref_name, "/database/cytoBand.txt.gz | \
gunzip -c | \
grep acen | \
cut -f-3 | \
sed 's/^chr//' | \
paste - - | \
cut -f1,2,6", sep=""), intern=T)
centromere <- read.table(textConnection(centromere_text), col.names=c("chr", "centr_begin", "centr_end"))

chrend_text <- system(paste("curl http://hgdownload.cse.ucsc.edu/goldenPath/", ref_name, "/database/chromInfo.txt.gz | \
gunzip -c | \
cut -f-2 | \
sed 's/^chr//'", sep=""), intern=T)
chrend <- read.table(textConnection(chrend_text), col.names=c("chr", "end"))

chr <- merge(chr, chrend, by.x="idx", by.y="chr")
chr <- merge(chr, centromere, by.x="idx", by.y="chr")


#VT
chr$name<-as.character(chr$name);
# this prints the NAs introduced by coercion message: X and Y cannot be replaced with a number
# commented by FGB
chr$idx<-as.numeric(as.character(chr$idx));
chr<-chr[order(chr$idx),c("end", "num_genes", "idx", "name", "centr_begin", "centr_end")];

gene$chr<-factor(as.character(gene$chr), 1:22, ordered=TRUE);
gene$gene<-as.character(gene$gene);
gene[gene$gene=="APOBEC3A_B", "gene"]<-"APOBEC3A.B" # old

gene[gene$gene=="GTF2H2C_2", "gene"]<-"GTF2H2C.2" #  added by wzk
gene[gene$gene=="C4B_2", "gene"]<-"C4B.2" #  added by wzk
gene[gene$gene=="TCONS_00029157", "gene"]<-"TCONS.00029157" #  added by wzk
gene[gene$gene=="XLOC_008559", "gene"]<-"XLOC.008559" #  added by wzk
gene[gene$gene=="XLOC_007697", "gene"]<-"XLOC.007697" #  added by wzk
gene[gene$gene=="XLOC_009911", "gene"]<-"XLOC.009911" #  added by wzk

# Are there any more genes with underscore in their name?
any(grepl('_', gene$gene))

# reindex chr
rownames(chr)<-NULL;
print(chr)

# order and reindex gene
gene <- gene[order(gene$gene),]
rownames(gene)<-NULL;

# split gene into chromosomes
chr_list<-gene[order(gene$chr, gene$begin, gene$end),]
chr_list<-split(data.frame(begin=chr_list$begin, end=chr_list$end, t=0, gene=chr_list$gene), chr_list$chr)

# save
system(paste("mkdir -p ./", ref_name, sep=""));
save(chr_list, gene, chr, file=paste("./", ref_name, "/data.RData", sep=""))
#now in your code you can do for example: hg18<-ref("human/hg18")
#end VT 

