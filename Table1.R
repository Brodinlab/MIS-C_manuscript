#Created by Giuseppe Rubens Pascucci (pascucci.1479790@studenti.uniroma1.it)
####################################################################

# ================================================== DATABASE ==================================================

mx <- read.table("DataIN/150720 clinical data.txt", sep = "\t", row.names = 1, header = T, check.names=F, stringsAsFactors = F)
mx$Gender <- gsub(" ", "", mx$Gender)
Group <- mx[,1]
Sex <- mx[,2]
mx <- as.matrix(mx[,-c(1,2)])
class(mx) <- "numeric"
tmx <- t(mx)



# ================================================== D'AGOSTINO-PEARSON FUNCTION ==================================================
dagostino.pearson.test <- function(x) {
  # from Zar (1999), implemented by Doug Scofield, scofield at bio.indiana.edu
  DNAME <- deparse(substitute(x))
  n <- length(x)
  x2 <- x * x
  x3 <- x * x2
  x4 <- x * x3
  # compute Z_g1
  k3 <- ((n*sum(x3)) - (3*sum(x)*sum(x2)) + (2*(sum(x)^3)/n)) /((n-1)*(n-2))
  g1 <- k3 / sqrt(var(x)^3)
  sqrtb1 <- ((n - 2)*g1) / sqrt(n*(n - 1))
  A <- sqrtb1 * sqrt(((n + 1)*(n + 3)) / (6*(n - 2)))
  B <- (3*(n*n + 27*n - 70)*(n+1)*(n+3)) / ((n-2)*(n+5)*(n+7)*(n+9))
  C <- sqrt(2*(B - 1)) - 1
  D <- sqrt(C)
  E <- 1 / sqrt(log(D))
  F <- A / sqrt(2/(C - 1))
  Zg1 <- E * log(F + sqrt(F*F + 1))
  # compute Z_g2
  G <- (24*n*(n-2)*(n-3)) / (((n+1)^2)*(n+3)*(n+5))
  k4 <- (((n*n*n + n*n)*sum(x4)) - (4*(n*n + n)*sum(x3)*sum(x)) - (3*(n*n - n)*sum(x2)^2) + (12*n*sum(x2)*sum(x)^2) - (6*sum(x)^4)) /(n*(n-1)*(n-2)*(n-3))
  g2 <- k4 / var(x)^2
  H <- ((n-2)*(n-3)*abs(g2)) / ((n+1)*(n-1)*sqrt(G))
  J <- ((6*(n*n - 5*n + 2)) / ((n+7)*(n+9))) * sqrt((6*(n+3)*(n+5)) /(n*(n-2)*(n-3)))
  K <- 6 + (8/J)*(2/J + sqrt(1 + 4/(J*J)))
  L <- (1 - 2/K) / (1 + H*sqrt(2/(K-4)))
  Zg2 <- (1 - 2/(9*K) - (L^(1/3))) / (sqrt(2/(9*K)))
  K2 <- Zg1*Zg1 + Zg2*Zg2
  pk2 <- pchisq(K2, 2, lower.tail=FALSE)
  RVAL <- list(statistic = c(K2 = K2), p.value = pk2, method = "D'Agostino-Pearson normality test\n\nK2 is distributed as Chi-squared with df=2", alternative = "distribution is not normal", data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}




# =========================================================== Differential Analysis ===========================================================

Group1 <- c("CoV2+","CoV2+", "CoV2+", "MIS-C", "MIS-C", "Kawasaki")
Group2 <- c("MIS-C", "Healthy", "Kawasaki", "Healthy", "Kawasaki", "Healthy")



t <- data.frame()

for(j in 1:length(Group1))
{
  
  x <- tmx[, Group %in% Group1[j]]
  y <- tmx[, Group %in% Group2[j]]
  
  smx <- rowSums(!is.na(x))
  smy <- rowSums(!is.na(y))
  
  ir <- which(smx < 3 | smy < 3)
  if(length(ir) > 0) {
    x <- x[-ir,]
    y <- y[-ir,]
  }
  
  a <- ncol(x)
  b <- ncol(y)
  m <- cbind(x, y)
  
  # Normality Test
  norX <- apply(m[, 1:a], 1, function(z){dagostino.pearson.test(z)$p.value})
  norY <- apply(m[, (a+1):(a + b)], 1, function(z){dagostino.pearson.test(z)$p.value})
  norX[is.na(norX)] <- 0
  norY[is.na(norY)] <- 0
  
  # Differential Test
  m <- cbind(m, norX, norY)
  pv <- apply(m, 1, function(z){ 
    if(z[a + b + 1] > 0.05 && z[a + b + 2] > 0.05){
      t.test(z[1:a], z[(a + 1):(a + b)], paired = F, na.action = "na.omit")$p.value
    } else { wilcox.test(z[1:a], z[(a + 1):(a + b)], paired = F, exact = F, na.action = "na.omit")$p.value }})
  adj <- p.adjust(pv, method="fdr")
  
  # Fold Change
  fc <- rowMeans(m[,1:a], na.rm = T) / rowMeans(m[,(a + 1):(a + b)], na.rm = T)
  
  # Final Table
  m <- cbind(Label = rownames(m), Group1 = Group1[j], Group2 = Group2[j], Comparison = paste0(Group1[j], " vs ", Group2[j]), pValue = pv, pAdj = adj, log2FoldChange = round(log2(fc + 0.01), 2))
  
  # if(class(m) %in% "matrix") t <- rbind(t, m)
  t <- rbind(t, m)
  
}

rm(x, y, a, b, m, j, norX, norY, pv, adj, fc)

write.table(t[,-c(2,3)], "DataOUT/Table1 Clinical Data (Diff All).txt", sep="\t", row.names = F, col.names = T,  quote = F)
t$pValue <- as.numeric(as.character(t$pValue))
t <- t[t$pValue < 0.05,]
write.table(t[,-c(2,3)], "DataOUT/Table1 Clinical Data (Diff DE).txt", sep="\t", row.names = F, col.names = T,  quote = F)




lab <- rownames(tmx)

smy <- cbind("Labels", "CoV2+", "MIS-C","Kawasaki", "Healthy")

for(j in 1:length(lab))
{
  
  a <- tmx[j, Group %in% "CoV2+"]
  b <- tmx[j, Group %in% "MIS-C"]
  c <- tmx[j, Group %in% "Kawasaki"]
  d <- tmx[j, Group %in% "Healthy"]
  
  la <- paste0(round(summary(a)[3],1), " (", round(summary(a)[2], 1), " - ", round(summary(a)[5], 1), ")")
  lb <- paste0(round(summary(b)[3],1), " (", round(summary(b)[2], 1), " - ", round(summary(b)[5], 1), ")")
  lc <- paste0(round(summary(c)[3],1), " (", round(summary(c)[2], 1), " - ", round(summary(c)[5], 1), ")")
  ld <- paste0(round(summary(d)[3],1), " (", round(summary(d)[2], 1), " - ", round(summary(d)[5], 1), ")")
  
  smy <- rbind(smy, cbind(lab[j], la, lb, lc, ld))
  
  
}


a <- table(Sex[Group %in% "CoV2+"])
b <- table(Sex[Group %in% "MIS-C"])
c <- table(Sex[Group %in% "Kawasaki"])
d <- table(Sex[Group %in% "Healthy"])
smy <- rbind(smy, cbind("Male : Female ", paste0(a[2], " : ", a[1]),  paste0(b[2], " : ", b[1]),  paste0(c[2], " : ", c[1]),  paste0(d[2], " : ", d[1])))

write.table(smy, "DataOUT/Table1 Clinical Data.txt", sep="\t", row.names = F, col.names = T,  quote = F)


tab <- rbind(cbind(a[1], a[2]), cbind(b[1], b[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

tab <- rbind(cbind(a[1], a[2]), cbind(c[1], c[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

tab <- rbind(cbind(a[1], a[2]), cbind(d[1], d[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

tab <- rbind(cbind(b[1], b[2]), cbind(c[1], c[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

tab <- rbind(cbind(b[1], b[2]), cbind(d[1], d[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)

tab <- rbind(cbind(c[1], c[2]), cbind(d[1], d[2]))
fisher.test(as.matrix(tab), y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)
