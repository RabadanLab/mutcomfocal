smooth<-function(mcf, w) UseMethod("smooth");
smooth.mcf_score<-function(mcf, w) with(ref(mcf), {
    y3<-colnames(mcf);
    x851 <- gene[mcf[, 1], 2:4];
    x851 <- data.frame(x851[, 1], (x851[, 2]+x851[, 3])/2, 1:nrow(x851));
    x852 <- order(x851[,1], x851[,2]);
    x851 <- x851[x852, ];
    r <- cbind(x851, mcf[x852, 2], rep(0, nrow(x851)), mcf[x852, 1]);
    
    if (w>0) {
        y2 <- (chr[, 5]+chr[, 6])/2;
        for(c in unique(r[,1])) {
            i<-as.numeric(c);
            for(j1 in 1:ifelse(is.na(y2[i]), 1, 2)) {
                l1 <- r[, 1]==c;
                if (!is.na(y2[i])) 
                    if (j1 == 1) {
                        ll <- l1 & r[, 2]<y2[i];
                    } else {
                        ll <- l1 & r[, 2]>=y2[i];
                    }
                
                x181 <- 1:sum(l1);
                x182 <- matrix(0, sum(l1),sum(l1));
                for(k in 1:sum(l1)) {
                    x182[,k] <- exp(-((x181-x181[k])^2)/w^2);
                    x182[,k] <- x182[,k]/sum(x182[,k]);
                }
                r[l1, 5] <- x182 %*% r[l1,4];
            }
        }
    } else {
        r[, 5] <- r[, 4];
    }

    r <- r[order(r[,3]), c(6,4:5)];
    colnames(r)<-c(y3, paste(y3[2], "smooth", sep="_"));
    class(r)<-c("mcf_score", class(r));
    ref(r)<-ref(mcf);
    r;
});
