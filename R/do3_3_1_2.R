do12<-function (x, y) {
    y[which.min(x)];
}


do3_3_1_2<-function() with(parent.frame(), {
    x8410 <- data[, 1]==0;
    x8411 <- data[, 1]==1;
    x8412 <- data[, 1]==2;

    x8416 <- x8132_1[x82, 1]>=0.3;
    x8417 <- x8132_1[x82, 1];

    q1 <- with(data, length(unique(data[data[, 1]<=1, 3])));
    q2 <- with(data, length(unique(data[data[, 1]==2, 3])));
    
    x8130[, 19] <- tapply(as.numeric(data[x8410, 3])*x8416[x8410], data[x8410, 4], function (x) length(unique(x[x!=0])))/q1;

    z1 <- tapply(x82[x8410], data[x8410, 4], function (x) do12((x8132_1[x, 5]*x8132_1[x, 3])/x8132_1[x, 1], x));
    z2 <- rbind(rep(NA, 3), x8132_1[, c(4,1,5)]);
    x8130[,23:25] <- z2[z1+1,];

    z1 <- tapply(x82[x8410], data[x8410, 4], function (x) do12(x8132_1[x, 6]/x8132_1[x, 1], x));
    z2 <- rbind(rep(NA,2), x8132_1[, c(1,6)]);
    x8130[, 26:27] <- z2[z1+1, ];

    x8130[, 28] <- tapply(x8132_1[x82[x8412], 7], data[x8412, 4], min);
    x8130[, 29] <- tapply(x8417[x8410],  data[x8410, 4], sum);

    x8130[is.na(x8130)] <- 0;

    x8131[, 19] <- tapply(as.numeric(data[x8411, 3])*x8416[x8411], data[x8411, 4], function (x) length(unique(x[x!=0])))/q1;

    z1 <- tapply(x82[x8411], data[x8411, 4], function (x) do12((x8132_1[x, 5]*x8132_1[x, 3])/x8132_1[x, 1], x));
    z2 <- rbind(rep(NA, 3), x8132_1[, c(4,1,5)]);
    x8131[, 23:25] <- z2[z1+1, ];
    z1 <- tapply(x82[x8411], data[x8411, 4], function (x) do12(x8132_1[x, 6]/x8132_1[x, 1], x));
    z2 <- rbind(rep(NA,2), x8132_1[, c(1,6)]);
    x8131[, 26:27] <- z2[z1+1, ];

    x8131[, 28] <- tapply(x8132_1[x82[x8412], 7], data[x8412, 4], min);
    x8131[, 29] <- tapply(x8417[x8411], data[x8411, 4], sum);

    x8131[is.na(x8131)] <- 0;

    x8133[, 2] <- tapply(data[x8412, 3], data[x8412, 4], function (x) length(unique(x[x!=0])))/q2;
    x8133[is.na(x8133[,2]), 2] <- 0;
    x8133[, 5] <- tapply(x8132_1[x82[x8412], 7], data[x8412, 4], min);
    x8133[is.na(x8133[,5]), 5] <- 0;    
})
