#Described below is the 'weighted area between the curves' (wABC; Edelen, Stucky, and Chandra, 2013)
#The code is applicable to 2-PL, 3-PL, and GRM models
#The example creates the expected scores curves and wABC from Figure 1 in Stucky et al. (2014)
#Comments highlight user inputs

theta <- seq(-5,5,.1)
RefPDF <- dnorm(theta, mean = 0, sd = 1, log = FALSE) #Reference group mean and SD 
RefPDF <- normalPDF/sum(normalPDF)
FocPDF <-  dnorm(theta, mean = -0.29, sd = 1.06, log = FALSE) #Focal group mean and SD
FocPDF <- FocPDF/sum(FocPDF)

Nref <- 421  #Reference group sample size
Nfoc <- 1825 #Focal group sample size

aj <- c(.53,1.18) #Enter the slope for reference group and the slope for focal group (there is no scaling constant in the model)
b1 <- c(-5.43,-3.25,-1.1,0.87)   #Enter the thresholds for reference group
b2 <- c(-2.45,-1.36,-0.19,0.91)   #Enter the thresholds for focal group
c <- c(0,0) #Enter guessing parameter for 3-PL

Ncats <- 5 #Enter the number of response categories for the item

bj <- list(b1,b2)
nitems <- length(aj)

Tlist <- NULL
for (i in 1:nitems) {
  tempi <- NULL
  for (j in 1:(length(bj[[i]])+1)) {
    tempj <- rep(0,length(theta))
    tempi <- c(tempi,list(tempj))
  }
  Tlist <- c(Tlist,list(tempi))
}
TCC <- Tlist

for (i in 1:nitems) {
  Tlist[[i]][[1]] <- (1- (c[i]+(1-c[i])/(1+exp(-aj[i]*(theta-bj[[i]][[1]])))))
  TCC[[i]] <- Tlist[[i]][[1]]*0
  Tlist[[i]][[(length(bj[[i]])+1)]] <- c[i]+(1-c[i])/(1+exp(-aj[i]*
                                                  (theta-bj[[i]][[length(bj[[i]])]])))
  TCC[[i]] <- TCC[[i]] + (Tlist[[i]][[(length(bj[[i]])+1)]])*(length(bj[[i]]))
  if  (Ncats > 2) {	
    for (j in 2:length(bj[[i]])) {
      Tlist[[i]][[j]] <- ((1/(1+exp(-aj[i]*(theta-bj[[i]][[j-1]]))))-
                            (1/(1+exp(-aj[i]*(theta-bj[[i]][[j]])))))
       TCC[[i]] <- TCC[[i]] + (Tlist[[i]][[j]])*(j-1)
    }
  }
}
wABC_ref <- sum(abs(TCC[[1]]-TCC[[2]])*RefPDF)
wABC_foc <- sum(abs(TCC[[1]]-TCC[[2]])*FocPDF)
wABC <- (wABC_foc*(Nfoc/(Nfoc+Nref)))+ (wABC_ref*(Nref/(Nfoc+Nref)))

plot(theta,TCC[[1]],type="l", xlim=c(-3,3), ylim=c(0,length(b1)), 
     xlab=expression(theta),ylab="Expected score", main="Expected score curves")
lines(theta,TCC[[2]],lty="dashed")
round(wABC, digits=2)
legend ((.6*min(theta)),1*length(b1), 
        c("Reference group","Focal group"), lty=c(1,2), lwd=c(2,2), 
        col=c("black","black"),bty="n",cex=c(1,1))
text(.4*min(theta),.8*length(b1), paste0 ("wABC =", wABC=(round(wABC,digits=2))))

#Edelen, M. O., Stucky, B. D., & Chandra, A. (2013). Quantifying ‘problematic’ DIF 
#within an IRT Framework: Application to a cancer stigma index. Quality of Life 
#Research, 1-9.

#Stucky, B. D., Edelen, M. O., Tucker, J. S., Shadel, W. G., Cerully, J., Kuhfeld, M., Hansen, M., 
#and Cai, L. (2014). Development of the PROMIS® Negative Psychosocial 
#Expectancies of Smoking item banks. Nicotine and Tobacco Research, 16(Suppl 3), 
#S232-S240.

