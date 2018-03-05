do3_3_3_4<-function() with(parent.frame(), {
    x82<-as.numeric(factor(do.call(paste, data[,c(1,3,5)])));
    x82<-match(x82, unique(x82))
    x83_1<-which(!duplicated(x82));
    x81<-data[x83_1, c(1,3,5)]
    x83<-matrix(NA, nrow(x81), 7);
    x83[,1]<-data[x83_1, 2]
    x83[,4]<-tapply(rep(1, length(x82)), x82, sum);

    x83_2<-as.numeric(factor(do.call(paste, x81[,1:2])));
    x83_2<-match(x83_2, unique(x83_2));
    x83_1<-which(!duplicated(x83_2));
    x81_1<-x81[x83_1, 1:2];
    x83_3<-tapply(x83[,1], x83_2, sum);
    x83_5<-tapply(rep(1, length(x83_2)), x83_2, sum);
    x83[,5]<-x83_3[x83_2];
    x83[,7]<-x83_5[x83_2];
})


