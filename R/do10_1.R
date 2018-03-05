do10_1<-function(x49_1) {
    z1 <- x49_1;
    z2 <- 1:length(x49_1);

    z2 <- z2[z1!=0];
    z1 <- z1[z1!=0];
    z1 <- z1/sum(z1);
    p <- 1; j <- 1;
    x49_2 <- rep(0, length(x49_1));
    x49_3 <- c();
    while (TRUE) {
        h1 <- exp(sum(z1*log(z1)));
        if ((sum(z1<h1)==length(z1)) && length(z1)>0)
            h1 <- min(z1);
        x49_2[z2[z1>=h1]] <- j;
        if (sum(z1<h1) == length(z1)) 
            break
        x49_3 <- rbind(x49_3, c(sum(z1[z1>=h1]), p, h1, sum(z1>=h1)));
        z2 <- z2[z1<h1];
        z1 <- z1[z1<h1]; 
        p <- p*sum(z1);
        z1 <- z1/sum(z1);
        j <- j + 1;
    }

    if (sum(x49_2==0)>0) {
        x49_3 <- rbind(x49_3, c(1, p, 0, sum(x49_2==0)));
        x49_2[x49_2==0] <- j;
    }

    list(x49_2=x49_2, x49_3=x49_3);
}
