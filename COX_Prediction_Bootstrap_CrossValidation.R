library("survival")
library("ggplot2")
library("ggpubr")
library("magrittr")
library("survminer")
library("splines")
library("dplyr")
library("BBmisc")
library("boot")

maxbootstraps <- 1000
maxiterations <- 1000

BM_t <-
  read.table(
    "/Users/mercedehmovassagh/Desktop/CorEBV7.txt",
    header = T,
    sep = "\t"
  )

runeverything <- function(train, test) {
  mel.cox <-
    coxph(
      Surv(time, status) ~ Age + sex + Eosinophil + Lymphocyte +	Myelocite + Neutrophil +
        Monocyte +M0+M1+M2+M3+M4+M5+M6+M7+
        trans8_21 + inv_16 + trans_8 + q_5d + q_7d + trans9_22 +
        trans4_11 + trans9_11 + trans15_17 + idh1_r140 +idh1_r172+
        idh1_r132 + NPMc + strata(EBV),
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
        Surv(time, status) ~ Age + sex + Eosinophil + Lymphocyte +	Myelocite + Neutrophil +
          Monocyte +M0+M1+M2+M3+M4+M5+M6+M7+idh1_r172+
          trans8_21 + inv_16 + trans_8 + q_5d + q_7d + trans9_22 +
          trans4_11 + trans9_11 + trans15_17 + idh1_r140 +
          idh1_r132 + NPMc + strata(EBV),
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


xval <- 1:nrow(BM_t)
ret1 <-
  c() # this will in the end contain nrow(BM_t) number of counts
ret2 <-
  c() # this will in the end contain nrow(BM_t) number of medians
#ret3 <-
  #c() # this will in the end contain nrow(BM_t) number of means
for (k in 1:nrow(BM_t)) {
  cat(sprintf("cross-validation round %d\n", k))
  train <-
    BM_t[xval != k, , drop = F] # training set is always nrow(BM_t) - 1 items
  test <- BM_t[xval == k, , drop = F] # test set is always 1 item
  preds  <- runeverything(train, test)
  # length(which(log(bla$t[, 1]) > 0))  ### counting how many log(pred) out of 1000 are > 0
  ret1 <- c(ret1, length(which(log(preds$t[, 1]) > 0)))
  ret2 <- c(ret2, median(log(preds$t[, 1])))
  #ret3 <- c(ret3, mean(log(preds$t[, 1])))
  print(ret1)
  print(ret2)
  #print(ret3)
}
# here we attach the counts, medians and means to BM_t (we know the order is the same!)
BM_t$poslogriskcount <- ret1
BM_t$logriskmedian <- ret2
BM_t$logriskmean <- ret3

E <- (BM_t$status)
P <- (BM_t$EBV)
K <- cbind(E, P)
write.table(E, "/Users/mercedehmovassagh/Desktop/trainSurv1", sep = "\t")
library(ROCR)

perf <-
  read.table(
    "/Users/mercedehmovassagh/Desktop/tryROCEBV.txt",
    header = T,
    sep = "\t"
  )

pred1 <- prediction(perf$Predictions, perf$labels)
perf1 <- performance(pred1, "tpr", "fpr")
plot(perf1,col = "blue",print.auc=TRUE)
abline(a = 0, b = 1)
text(x = 0.8, y = 0.2, paste("AUC=0.72"), 
      cex = 1, col = "black")

auc.perf = performance(pred1, measure = "auc")
auc.perf@y.values


#######Median
perf <-
  read.table(
    "/Users/mercedehmovassagh/Desktop/tryROCEBV.txt",
    header = T,
    sep = "\t"
  )

pred1 <- prediction(perf$Predictions, perf$labels)
perf1 <- performance(pred1, "tpr", "fpr")
plot(perf1,col = "blue",print.auc=TRUE)
abline(a = 0, b = 1)
text(x = 0.8, y = 0.2, paste("AUC=0.74"), 
     cex = 1, col = "black")

auc.perf = performance(pred1, measure = "auc")
auc.perf@y.values
