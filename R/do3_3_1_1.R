do3_3_1_1<-function() with(parent.frame(), {
    x9 <- c(0.5, 0.5);

    if (!exists("x9_1_1")) x9_1_1 <- matrix(1, nrow(ref(data)$gene), 3);
    if (!exists("x9_1_2")) x9_1_2 <- matrix(1, nrow(ref(data)$gene), 3);

    x9_1_1 <- sweep(x9_1_1, 2, colSums(x9_1_1), FUN="/");
    x9_1_2 <- sweep(x9_1_2, 2, colSums(x9_1_2), FUN="/");

    x9_2 <- cbind(c(x9_1_1)[as.numeric(data[, 4])+nrow(x9_1_1)*data[, 1]], 
                  c(x9_1_2)[as.numeric(data[, 4])+nrow(x9_1_2)*data[, 1]]);
    
    x83[, 2] <- tapply(x9_2[,1], x82, sum);
    x83[, 3] <- tapply(x9_2[,2], x82, sum);

    x83_4 <- tapply(x83[,1]*x83[,2], x83_2, sum);
    x83[, 6] <- x83_4[x83_2];

    x8410 <- data[, 1]==0;
    x8411 <- data[, 1]==1;
    x8412 <- data[, 1]==2;
    
    x842 <- 1/x83[x82, 4];
    
    x843 <- x83[x82, 1]/x83[x82, 4];

    x844 <- x83[x82, 1]/(x83[x82, 7]*x83[x82, 4]);

    x8414 <- (x9_2[, 1]*x83[x82, 1])/x83[x82, 6];
    x84140 <- x8414;
    x84140[x8412] <- x84140[x8412]*x9[1];
    x84141 <- x8414;
    x84141[x8412] <- x84141[x8412]*x9[2];
    
    x8415 <- (x9_2[, 2]*x83[x82, 1])/(x83[x82, 5]*x83[x82, 3]);
    x84150 <- x8415;
    x84150[x8412] <- x84150[x8412]*x9[1];
    x84151 <- x8415;
    x84151[x8412] <- x84151[x8412]*x9[2];

    x8132_1 <- x83;

    x8130<-matrix(NA, length(levels(data[,4])), 29);
    x8130[, 1] <- as.numeric(levels(data[,4]));
    x8130[, 2] <- tapply(x842[x8410], data[x8410, 4], sum);
    x8130[, 3] <- tapply(x843[x8410], data[x8410, 4], sum);
    x8130[, 4] <- tapply(x844[x8410], data[x8410, 4], sum);
    x8130[, 5] <- tapply(x8414[x8410], data[x8410, 4], sum);
    x8130[, 6] <- tapply(x8415[x8410], data[x8410, 4], sum);
    x8130[, 7] <- tapply(x84140[!x8411], data[!x8411, 4], sum);
    x8130[, 8] <- tapply(x84150[!x8411], data[!x8411, 4], sum);
    x8130[, 9] <- tapply(x8414[x8412], data[x8412, 4], sum);
    x8130[is.na(x8130)] <- 0;
    x8130[, 2:9] <- sweep(x8130[, 2:9], 2, colSums(x8130[,2:9]), FUN="/");
    x8130[, 13] <- x8130[, 5]*x8130[, 6]; x8130[, 13] <- x8130[, 13]/sum(x8130[, 13]);
    
    x8131<-matrix(NA, length(levels(data[,4])), 29);
    x8131[, 1] <- as.numeric(levels(data[,4]));
    x8131[, 2] <- tapply(x842[x8411], data[x8411, 4], sum);
    x8131[, 3] <- tapply(x843[x8411], data[x8411, 4], sum);
    x8131[, 4] <- tapply(x844[x8411], data[x8411, 4], sum);
    x8131[, 5] <- tapply(x8414[x8411], data[x8411, 4], sum);
    x8131[, 6] <- tapply(x8415[x8411], data[x8411, 4], sum);
    x8131[, 7] <- tapply(x84141[!x8410], data[!x8410, 4], sum);
    x8131[, 8] <- tapply(x84151[!x8410], data[!x8410, 4], sum);
    x8131[, 9] <- tapply(x8414[x8412], data[x8412, 4], sum);
    x8131[is.na(x8131)] <- 0;
    x8131[, 2:9] <- sweep(x8131[, 2:9], 2, colSums(x8131[,2:9]), FUN="/");
    x8131[, 13] <- x8131[, 5]*x8131[, 6]; x8131[, 13] <- x8131[, 13]/sum(x8131[, 13]);
    
    x8133<-matrix(NA, length(levels(data[,4])), 6)
    x8133[, 1] <- tapply(x8414[x8412], data[x8412, 4], sum);
    x8133[is.na(x8133)] <- 0;
    x8133[, 1] <- x8133[, 1]/sum(x8133[, 1]);

    x8130[, 11] <- 0.5*(x8130[, 5]+x8130[, 6])*x8133[, 1]; x8130[, 11] <- x8130[, 11]/sum(x8130[, 11]);
    x8131[, 11] <- 0.5*(x8131[, 5]+x8131[, 6])*x8133[, 1]; x8131[, 11] <- x8131[, 11]/sum(x8131[, 11]);
    
    x8133[, 6] <- 0.25*(x8130[, 11]+x8130[, 13]+x8131[, 11]+x8131[, 13]);
})
