"[.mcf_score"<-function(x, ...) { r<-NextMethod("[", x); ref(r)<-ref(x); r }
"[<-.mcf_score"<-function(x, i, j, value) { r<-NextMethod("[<-", x); ref(r)<-ref(value); r }

as.score<-function(x, f) UseMethod("as.score", x);
as.score.data.frame<-function(x, f) {
    ref(x)<-f;
    class(x)<-c("mcf_score", "data.frame");
    x;
}

score<-function(data, perm=0, iter=1, prior_recurr=NULL, prior_focal=NULL) UseMethod("score");
score.mcf_lesions<-function(data, perm=0, iter=1, prior_recurr=NULL, prior_focal=NULL) {
    options(stringsAsFactors=FALSE);
    if (is.null(prior_recurr)) prior_recurr <- matrix(1, nrow(ref(data)$gene), 3);
    if (is.null(prior_focal)) prior_focal <- matrix(1, nrow(ref(data)$gene), 3);
    data[,4]<-factor(data[,4]);
    x9_3<-iter;
    if (perm>0) {
        x9_4<-perm;
        x1_1 <- data;
        z2 <- levels(data[, 4]);
        x600 <- list();
        for(i1 in 1:x9_4) {
            cat(i1, "\n");

            for(i2 in unique(data[, 3])) {
                z1 <- z2[sample(1:length(z2), length(z2), replace=FALSE)];
                data[data[, 3] == i2, 4] <- z1[data[data[, 3] == i2, 4]];
            }
            
            x9_1_1 <- prior_recurr;
            x9_1_2 <- prior_focal;
            do3_3_3_4();
            do3_3_3_3();

            x8130<-x8130[,c(5,6,11,13)];
            colnames(x8130)<-paste("amp", c("recurr", "focal", "mut", "comb"), sep="_");
            x8131<-x8131[,c(5,6,11,13)];
            colnames(x8131)<-paste("del", c("recurr", "focal", "mut", "comb"), sep="_");
            x8133<-x8133[,c(1,6)];
            colnames(x8133)<-c("mut", "all");
            x8134<-cbind(x8130, x8131, x8133);
        
            x600[[i1]] <- x8134;
            
            x1 <- x1_1;
        }
        
        x8135<-matrix(NA, length(levels(data[,4])), 10);
        x601<-matrix(NA, length(levels(data[,4])), length(x600));
        for(i2 in 1:10) {
            for(i1 in 1:length(x600))
                x601[,i1]<-x600[[i1]][,i2];
            x8135[,i2]<-rowMeans(x601);
        }
        colnames(x8135)<-c(paste("amp", c("recurr", "focal", "mut", "comb"), sep="_"),
                            paste("del", c("recurr", "focal", "mut", "comb"), sep="_"),
                            c("mut", "all"));

        x8135<-as.data.frame(x8135);
        class(x8135)<-c("mcf_perms", class(x8135));
        x8135;
    } else {
        do3_3_3_4();
        do3_3_3_3();
        x8130 <- x602[[1]][[1]]; x8131 <- x602[[1]][[2]]; x8133 <- x602[[1]][[3]];
        do3_3_1_2();
#        x8130<-x8130[,c(1,5,6,11,13,19,23,24,25,26,27)];
#        colnames(x8130)<-c("idx", paste("amp", c("recurr", "focal", "mut", "comb", "freq",
#                                                 "min1_genes", "min1_cnv", "min1_cnv_nbrs", "min2_cnv", "min2_tot"), sep="_"));
#        x8131<-x8131[,c(5,6,11,13,19,23,24,25,26,27)];    
#        colnames(x8131)<-paste("del", c("recurr", "focal", "mut", "comb", "freq",
#                                        "min1_genes", "min1_cnv", "min1_cnv_nbrs", "min2_cnv", "min2_tot"), sep="_");
        x8130<-x8130[,c(1,5,6,11,13)];
        colnames(x8130)<-c("idx", paste("amp", c("recurr", "focal", "mut", "comb"), sep="_"));
        x8131<-x8131[,c(5,6,11,13)];
        colnames(x8131)<-paste("del", c("recurr", "focal", "mut", "comb"), sep="_");
        
#        x8133<-x8133[,c(1,2,5,6)];
#        colnames(x8133)<-c("mut", "mut_freq", "mut_min_nbrs", "all");
        x8133<-x8133[,c(1,6)];
        colnames(x8133)<-c("mut", "all");        
        x8134<-cbind(x8130, x8131, x8133);
        x8134<-as.data.frame(x8134);
        class(x8134)<-c("mcf_score", class(x8134));
        ref(x8134)<-ref(data);
        x8134;
    }

}
