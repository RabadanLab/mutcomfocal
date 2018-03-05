segsmerge<-function (l, fl_o=NULL, fl_i=FALSE) {
    z1<-1:2;
    cols<-c("b", "e");
    if (fl_i) {
        z1<-c(z1, 3);
        cols<-c("d", cols);
    }
    if (!is.null(fl_o)) {
        c<-0;
        z1<-c(z1, max(z1)+1);
        c<-c(cols, "t");
    }
    z2<-lapply(1:length(l), function (i) {
        r<-data.frame(l[[i]][,z1], i, 1:nrow(l[[i]]));
        colnames(r)<-c(cols, "i", "n");
        r;        
    });
    z2<-do.call(rbind, z2);
    if (fl_i) {
        z2<-z2[order(z2$d, z2$b, z2$e),];
    } else {
        z2<-z2[order(z2$b, z2$e),];
    }

    i<-1;
    r<-z2[i,,drop=FALSE];
    fl<-TRUE;
    fl2<-TRUE;
    q <- list();
    z5<-list();
    ll<-0;
    while (TRUE) {
        if (i<=nrow(z2) && fl2) {
            while (TRUE) {			
                if (length(q)==0) {
                    if (fl_i) d <- r$d;
                    b <- r$b;
                    e <- r$e;
                } else if (r$b == b) { 
                    e <- min(e, r$e);
                } else if (fl) {
                    fl <- FALSE;
                    e <- min(e, r$b-1);
                    break;
                } else {
                    ne <- r$e;
                    j<-1;
                    while (j<=length(q)) {
                        if (q[[j]]$e == e) {
                            if (!is.null(fl_o)) c <- c-(q[[j]]$t == fl_o);
                            q[[j]]<-NULL;
                        } else {
                            ne <- min(ne, q[[j]]$e);
                            j <- j+1;
                        }
                    }
                    if (r$b > e+1 && length(q)>0) {
                        b <- e + 1;
                        e <- min(ne, r$b-1); 
                        break; 
                    } else {
                        fl <- TRUE;
                        b <- r$b;
                        e <- ne;
                    }
                }
				
                if (!is.null(fl_o)) c <- c + (r$t == fl_o);
                q[[length(q)+1]]<-r;
                    
                i <- i+1;
                if (i<=nrow(z2)) {
                    r<-z2[i,,drop=FALSE];
                    fl2 = !fl_i || fl_i&&(r$d==d);
                }
                if (i>nrow(z2) || !fl2) break;
            }
        } else {
            fl <- TRUE;
            b <- e + 1;
            j<-1;
            while (j<=length(q)) 
                if (q[[j]]$e == e) {
                    if (!is.null(fl_o)) c <- c-(q[[j]]$t == fl_o);
                    q[[j]]<-NULL;
                } else { 
                    if (fl) {
                        ne <- q[[j]]$e;
                        fl <- FALSE;
                    } else {
                        ne <- min(ne, q[[j]]$e);
                    }
                    j <- j+1;
                }
            if (length(q)==0) {
                if (i>nrow(z2)) {
                    break;
                } else {
                    fl2 <- TRUE;
                }
            } else {
                e <- ne;
            }
        }
            
        if (length(q)>0 && (is.null(fl_o) || !is.null(fl_o)&&c>0)) {
            z3<-list(b=b, e=e, n=length(q));
            if (fl_i) z3<-c(d=d, z3);
            qd<-do.call(rbind, lapply(q, function (x) c(x$i, x$n)));
            colnames(qd)<-c("i", "n");
            ll<-ll+1;
            z5[[ll]]<-list(z3, qd);
        }
    }
    z5;
}

lesions<-function(r, cnv, mut) UseMethod("lesions");
lesions.mcf_ref<-function(r, cnv, mut) {
    options(stringsAsFactors=FALSE);
    require(hash);
    require(data.table);

    x2<-sort(unique(c(cnv[,1], mut[,1])));
    x2<-data.frame(x2, 1:length(x2));
    x3<-data.frame(cnv[,-1], cnv[,1]);
    x3[,1]<-gsub("24", "Y", gsub("23", "X", x3[,1]));
    x3<-data.frame(x3[order(x3[,1], x3[,2], x3[,3]),], 1:nrow(x3));

    chr<-list();
    for(c in sort(unique(x3[,1]))) {
        chr[[c]]<-list();
        chr[[c]]$x1<-data.frame(0, x3[x3[,1]==c,-1]);
    }
    
    g.x9.2<-data.frame(mut[,-1], mut[,1], 1:nrow(mut));
    x9.3<-sort(unique(g.x9.2[,1]));
        
    for(c in names(chr)) chr[[c]]<-within(chr[[c]], {
        cat(c, "\n");
        x9.7<-r$chr_list[[c]];
        x9.1<-data.frame(paste(x9.7[,4], c, sep="@"), x9.7[,4], x9.7[,1], x9.7[,2]);
        x9.2<-sort(unique(x9.1[,1]));
        
        x9.3<-g.x9.2;
        x9.3<-x9.3[x9.3[,1] %in% x9.2,];
        x9.4<-data.frame(x9.3, x9.1[match(x9.3[,1], x9.1[,1]),-1]);
        x9.4<-data.frame(1, x9.4[,4], x9.4, 1, x9.4[,3]);    
        x9.4<-x9.4[,c(1,2,4,7:ncol(x9.4))];
        x9.5<-data.frame(x9.1[,1], c, x9.1[,2:3]);
        
        x2<-data.frame(x1[,2:3], 1, x1[,7], 0+(x1[,5]<0), abs(x1[,5]), x1[,6]);
        
        x2.2.1<-data.frame(paste(c, 0, x2[,5], sep="_"), c, x2[,2:3], x2[6:7]);
        colnames(x2.2.1)<-1:6;
        x2.2.2<-data.frame(paste(c, x9.4[,1], x9.4[,7], sep="_"), c, x9.4[,4:5], 1+x9.4[,1], x9.4[,6]);
        colnames(x2.2.2)<-1:6;    
        x2.2<-rbind(x2.2.1, x2.2.2);
        
        x3.1<-segsmerge(list(x2, x9.7));
        
        f<-hash();
        x2_9<-as.matrix(x2[,c(1:2,5)]);
        x2_10<-as.matrix(x9.7[,1:2]);
        x3<-lapply(x3.1, function (x) {
            x2_1<-x[[2]][x[[2]][,"i"]==1,"n"];
            x2_2<-x[[2]][x[[2]][,"i"]==2,"n"];
            if (length(x2_1)>0 && length(x2_2)>0) {
                x3_3<-as.matrix(expand.grid(x2_1, x2_2));
                x3_3<-unique(x3_3);
                x3_4<-paste(x3_3[,1], x3_3[,2]);
                x3_5<-!has.key(x3_4, f);
                x3_4<-x3_4[x3_5];
                x3_3<-x3_3[x3_5,,drop=FALSE];
                x3_3<-cbind(x3_3, x2_9[x3_3[,1],,drop=FALSE], x2_10[x3_3[,2],,drop=FALSE]);
                x3_5<-x3_3[,5]==1 | x3_3[,5]!=1&x3_3[,3]<=x3_3[,6]&x3_3[,7]<=x3_3[,4];
                x3_4<-x3_4[x3_5];
                if (length(x3_4)>0) f[x3_4]<-1;
                as.data.frame(x3_3[x3_5,1:2,drop=FALSE]);
            }
        });
        x3<-as.matrix(rbindlist(x3));
        x3<-data.frame(x9.7[x3[,2], c(4,1,2)], paste(c, 0, x2[x3[,1],4], sep="_"), x2[x3[,1],5:7]);
        setnames(x3, as.character(1:ncol(x3)));
        
        x3.2<-data.frame(x9.4[,c(2,4,5)], paste(c, x9.4[,1], x9.4[,7], sep="_"), x9.4[,1]+1, x9.4[,c(6,3)]);
        colnames(x3.2)<-1:ncol(x3.2);
        
        x7.0<-rbind(x3, x3.2);
        x7.0_1<-paste(x7.0[,1], x7.0[,4], sep="@");
        x7.0_1<-do.call(rbind, strsplit(x7.0_1, "_", fixed=TRUE));
        x7.0<-data.frame(x7.0_1[,1], x7.0);
        x7.0<-data.frame(x7.0[,8], x7.0);
        x7.0<-data.frame(x7.0[,c(1,2,6,7,8)]);
    });
    x7.0<-do.call(rbind, lapply(chr, function (x) x$x7.0));

    x7.2<-data.frame(paste(r$gene[,1], r$gene[,2], sep="@"), 1:nrow(r$gene));
    
    x7.3<-sort(unique(x7.0[,3]))
    x7.3<-data.frame(x7.3, 1:length(x7.3));    
    x7.1<-data.frame(x7.0, x2[match(x7.0[,1], x2[,1]),-1])[,-1];
    x7.4<-data.frame(x7.1, x7.2[match(x7.1[,1], x7.2[,1]),2])[,-1];
    x7<-data.frame(x7.4, x7.3[match(x7.4[,1], x7.3[,1]),-1])[,-1];
    colnames(x7)<-1:ncol(x7);
    rownames(x7)<-1:nrow(x7);
    
    x7.5<-do.call(rbind, lapply(chr, function (x) x$x2.2));
    
    x8.1<-x2[,1];
    x7[,3]<-factor(x2[x7[,3],1], levels=x2[,1]);
    colnames(x7)<-c("type", "d", "sample", "gene", "id");
    class(x7)<-c("mcf_lesions", class(x7));
    ref(x7)<-r;
    x7
}
    
