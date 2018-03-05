do3_3_3_3<-function() with(parent.frame(), {
    x602 <- list();
    for(i3 in 1:x9_3) {
        cat(i3, "\n");
        do3_3_1_1();
        x9_1_1 <- cbind(x8130[, 5], x8131[, 5], x8133[, 1]);
        x9_1_2 <- cbind(x8130[, 6], x8131[, 6], x8133[, 1]);
        
        x8130[,-1] <- sweep(x8130[,-1], 2, colSums(x8130[,-1]), FUN="/");
        x8131[,-1] <- sweep(x8131[,-1], 2, colSums(x8131[,-1]), FUN="/");
        x8133[,-1] <- sweep(x8133[,-1], 2, colSums(x8133[,-1]), FUN="/");
        
        x602[[i3]] <- list(x8130, x8131, x8133);
    }
});
