library(ROCR)

#ROC <- c("~/ret_inv16_1.txt", "~/ret_Tri_1.txt", "~/ret_1.txt")

ROC <- c("~/ret_EBV_1.txt", "~/ret_woEBV_1.txt")

#names_legend<-c("Inversion 16, AUC=%.2f", "Trisomy 21, AUC=%.2f", "EBV, AUC=%.2f")

names_legend<-c("Model With EBV Addition, AUC=%.2f", "Model Without EBV, AUC=%.2f")

##If more add instead if 1:3-1:n
pdfname<-"Roc.pdf"
#Colors
#colorROC<-c("#9e6ebd","#7aa457","#cb6751")

colorROC<-c("#9e6ebd","#7aa457")


##add=True is doing multiple lines
pdf(pdfname, width = 7, height = 7)
AUCs<- c()
###Depending on how many lines u want 1:2,1:3,1:4
#for (i in 1:3) {
for (i in 1:2) {
  BM_t<-read.table(ROC[i],header=TRUE)
  pred1 <- prediction(BM_t$poslogriskcount, BM_t$status)
  perf1 <- performance(pred1, "tpr", "fpr")
  if (i==1) {
    plot(perf1, col = colorROC[i],lwd=3 , print.auc = TRUE, xlim=c(0,1), ylim=c(0,1))
  } else {
    plot(perf1, col = colorROC[i], lwd=3, print.auc = TRUE, add=T, xlim=c(0,1), ylim=c(0,1))
  }
  auc.perf = performance(pred1, measure = "auc")
  AUCs <- c(AUCs, auc.perf@y.values[[1]])
}

abline(a = 0, b = 1)
legend("bottomright",  sprintf(names_legend,AUCs), col=colorROC, lwd=5, cex=1.3)


dev.off()