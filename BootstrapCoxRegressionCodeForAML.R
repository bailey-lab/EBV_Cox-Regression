library("survival")
library("ggplot2")
library("ggpubr")
library("magrittr")
library("survminer")
library("splines")
library("dplyr")
library("BBmisc")
library("boot")
library("doParallel")
library("foreach")
library("ROCR")

registerDoParallel(cores=7)

maxbootstraps <- 1000
maxiterations <- 1000


for (ROUND in 1:1) {
  cat(sprintf("ROUND: %d\n", ROUND))
  
  BM_t <-
    read.table(
      "/Users/mercedehmovassagh/Desktop/CorEBV7.txt",
      header = T,
      sep = "\t"
    )
  
  pdfname <- sprintf("roc_%d.pdf", ROUND)
  retname <- sprintf("ret_%d.txt", ROUND)
  aucname <- sprintf("auc_%d.txt", ROUND)
  #If you wanted to check random one use bellow line if not just take it out like now!
  #BM_t$EBV <- sample(BM_t$EBV)
  
  runeverything <- function(train, test) {
    mel.cox <-
      coxph(
        Surv(time, status) ~ EBV+ Age + Eosinophil+ inv_16+Tri21,
        data = train,
        iter.max = maxiterations
      )
    
    agec <- cut(train$Age, c(0, 60, 100))
    
    mel.cens <-
      survfit(Surv(time - 0.001 * (status == 2), status != 2) ~ strata(agec), data = train)
    mel.surv <- survfit(mel.cox)
    mel.fun <- function(d) {
      cox <-
        coxph(
          Surv(time, status) ~ EBV+ Age + Eosinophil+ inv_16+Tri21,
          data = d,
          iter.max = maxiterations
        )
      
      v <- predict(cox,
                   newdata = test,
                   type = "risk")
      return(v)
    }
    
    mel.str <- cbind(train$EBV, agec)
    mel.mod <- censboot(
      train,
      mel.fun,
      R = maxbootstraps,
      F.surv = mel.surv,
      G.surv = mel.cens,
      cox = mel.cox,
      strata = mel.str,
      sim = "ordinary"
    )
    return(mel.mod)
  }
  
  N <- nrow(BM_t)
  xval <- 1:N
  #ret1 <- c() # this will in the end contain nrow(BM_t) number of counts
  #ret2 <- c() # this will in the end contain nrow(BM_t) number of medians
  #ret3 <-
  #c() # this will in the end contain nrow(BM_t) number of means
  #for (k in 1:N)
  #N <- 3
  #ret <- lapply(1:N, function(k) {
  ret <- foreach(k=1:N) %dopar% {
    cat(sprintf("cross-validation round %d\n", k))
    train <-
      BM_t[xval != k, , drop = F] # training set is always nrow(BM_t) - 1 items
    test <- BM_t[xval == k, , drop = F] # test set is always 1 item
    preds  <- runeverything(train, test)
    # length(which(log(bla$t[, 1]) > 0))  ### counting how many log(pred) out of 1000 are > 0
    ret1 <- length(which(log(preds$t[, 1]) > 0))
    ret2 <- median(log(preds$t[, 1]))
    #ret3 <- c(ret3, mean(log(preds$t[, 1])))
    #print(ret1)
    #print(ret2)
    #print(ret3)
    c(ret1, ret2)
  }
  ret <- do.call(rbind.data.frame, ret)
  names(ret) <- c("poslogriskcount", "logriskmedian")
  BM_t <- cbind(BM_t, ret)
  
  # make the files
  write.table(BM_t, quote=F, col.names=T, row.names=F, file=retname)
  
  
  
  
  
  ###Generate AUC for everything
  pred1 <- prediction(BM_t$poslogriskcount, BM_t$status)
  perf1 <- performance(pred1, "tpr", "fpr")
  pdf(pdfname, width=7, height=7)
  plot(perf1,col = "blue",print.auc=TRUE)
  abline(a = 0, b = 1)
  auc.perf = performance(pred1, measure = "auc")
  myauc <- auc.perf@y.values[[1]]
  text(x = 0.8, y = 0.2, sprintf("AUC=%.2f", myauc), 
       cex = 1, col = "black")
  dev.off() # close the pdf
  
  write.table(data.frame(auc=myauc), quote=F, col.names=T, row.names=F, file=aucname)
  
}