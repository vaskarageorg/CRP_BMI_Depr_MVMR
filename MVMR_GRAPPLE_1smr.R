#Load the GRAPPLE package
library(GRAPPLE)
grapple_coll = function(X1,X2,G1,G2,Y){
#Collider Correction
XGdata2 = data.frame(X1, G1,G2)             
FIT2 = summary(lm(XGdata2))             # estimated linear association
BetaX1G = FIT2$coef[-1,1]; seBetaX1G = FIT2$coef[-1,2] 
########################
YXGdata          = data.frame(Y, X1,X2,G1,G2)
FIT4             = summary(lm(YXGdata))                  # estimated associations
alphahatstar = FIT4$coef[-c(1,2,3),1]; se.alphahatstar  = FIT4$coef[-c(1,2,3),2]
betastarX1 = FIT4$coef[2,1]; se.betastarX1 = FIT4$coef[2,2] 
betastarX2 = FIT4$coef[3,1]; se.betastarX2 = FIT4$coef[3,2] 

#
X2Gdata2 = data.frame(X2, G1,G2)             
FIT2 = summary(lm(X2Gdata2))             # estimated linear association
BetaX2G = FIT2$coef[-1,1]; seBetaX2G = FIT2$coef[-1,2] 

data_GR1 = data.frame(gamma_exp1 = BetaX1G, gamma_exp2 = BetaX2G,
                      se_exp1 = seBetaX1G, se_exp2 = seBetaX2G,
                      gamma_out1 = alphahatstar, se_out1 = se.alphahatstar)
#MR GRAPPLE allows for an explicit correction when the exposures X1 and X2 are highly correlated.
#We therefore construct such a matrix for our purpose. 
cordf21=data.frame(X1, X2, Y); cordf21 = cordf21[complete.cases(cordf21),]
cor.mat = cor(cordf21); cor.mat[3,c(1,2)]=0; cor.mat[c(1,2),3]=0

Fit= grappleRobustEst(data = data_GR1,cor.mat = cor.mat,plot.it = F)
betaestX1 = betastarX1 + Fit$beta.hat[1]; sebetaestX1 = sqrt(se.betastarX1^2 + diag(Fit$beta.var)[1])
betaestX2 = betastarX2 + Fit$beta.hat[2]; sebetaestX2 = sqrt(se.betastarX2^2 + diag(Fit$beta.var)[2])

##get conditional instrument strength
#gencov2 = list(); Gall = data.frame(G1,G2); Xall = data.frame(X1,X2)
#for (j in 1:length(c(BetaX1G))){
#  gencov2[[j]]=matrix(nrow = ncol(Xall),ncol = ncol(Xall))
#  for (k in 1:ncol(Xall)){
#    resid_uk = (lm(Xall[,k]~Gall[,j])$resid)
#    for (m in 1:ncol(Xall)){
#      resid_um = (lm(Xall[,m]~Gall[,j])$resid)
#      gencov2[[j]][k,m] = (solve(((t(Gall[,j]) %*%Gall[,j])))/nrow(Gall)) * sum(resid_uk*resid_um)
#    }
#  } }

#CFS_1 = MVMR::strength_mvmr(MVMR::format_mvmr(BXGs = data.frame(BetaX1G,BetaX2G),
#                                                BYG = alphahatstar,
#                                                seBXGs = data.frame(seBetaX1G,seBetaX2G),
#                                                seBYG = se.alphahatstar),
#                            gencov = gencov2)
return(list(betaestX1=betaestX1,sebetaestX1=sebetaestX1,betaestX2=betaestX2,sebetaestX2=sebetaestX2 ))
}
