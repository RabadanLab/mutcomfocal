peaks<-function(mcf) UseMethod("peaks");
peaks.mcf_score<-function (mcf) with(ref(mcf), {
    options(stringsAsFactors=FALSE);
    x2<-data.frame(mcf, ref(mcf)$gene[mcf$idx,]);
    x2<-data.frame(x2[order(x2$chr, x2$begin),], 1:nrow(x2));

    res<-do10_2(x2[,2], x2$chr);
    x2<-cbind(x2, tier=res$x49_2);
    x2<-x2[x2[,8]>0,];
    x2<-cbind(x2, cumsum(c(TRUE, tail(x2$chr, -1)!=head(x2$chr, -1)) | c(FALSE, head(x2[,7],-1)+1!=tail(x2[,7], -1))));
    x2<-x2[order(x2[,9], -x2[,2]),];
    x3<-x2[!duplicated(x2[,9]),2];
    x2<-data.frame(x2, x3[x2[,9]]);
    x2<-x2[x2[,2]==x2[,10],];
    x5<-table(x2[,9]);
    x5<-ifelse(x5>1, paste("(+", x5-1, ")", sep=""), "");
    x6<-sapply(split(x2$gene, x2[,9]), function (x) paste(x, collapse=","));
    x7<-paste(paste(x2$gene, x2$chr, sep="@")[!duplicated(x2[,9])], x5, sep="");
    x8_1<-tapply(x2[,5], x2[,9], min);
    x8_2<-tapply(x2[,6], x2[,9], max);    
    x4<-data.frame(x2[!duplicated(x2[,9]),c(1:2,8,4)], begin=x8_1, end=x8_2, name=x7, genes=x6);
    x4<-x4[order(x4[,1]),];
    x4
})
