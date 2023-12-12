

# On peut ajouter ici des codes comun pour que ça soit plus facile





##### PARTIE GRAPHES Issu du code du prof


# ---------------------------------------------------
# Traceplot (see also ?traceplot - coda package)
# ---------------------------------------------------
varnames(mcmc)
windows()
par(mfrow=c(2,3))
traceplot(mcmc[,'K'],ylab="K")
traceplot(mcmc[,'r'],ylab="r")
traceplot(mcmc[,'q'],ylab="q")
traceplot(mcmc[,'C_MSY'],ylab="C_MSY")
traceplot(mcmc[,'h_MSY'],ylab="h_MSY")
traceplot(mcmc[,'B_MSY'],ylab="B_MSY")
traceplot(mcmc[,'sigma2p'],ylab="sigma2p")


# ----------------------------------------------
# Catches - predicted vs observed
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

par(mfrow = c(1,1), mar=c(5,5,1,1))

col = c("lightgray","gray")

# x = "C_pred"
x = "C"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
n_plot <- dim(X)[2]

q.mean <- apply(X,2,mean)
q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
  q.inf[k] <- quantile(X[,k],probs=0.05)
  q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
  q.inf1[k] <- quantile(X[,k],probs=0.25)
  q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=Year[1:n_plot],y=q.sup*1.2,ylim=c(0,max(q.sup*1.2)),ylab="Catches",type="n",las=1, main = "Fit Catches")
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup,rev(q.inf)),col = col[1],border=NA)
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=Year[1:n_plot], y=q.mean, col="red", lwd=2, type="l")
points(Year[1:n_plot],data$C_obs[1:n_plot],pch=20,cex=1,col="blue",type="p")

legend(	legend = c("Cred interv. (95%)","C pred.","C obs."),
        col = c("lightgray","red","blue"), 	
        lty = c(3,2,NA), pch=c(NA,NA,20), lwd = c(2,2,NA), x = "topright", cex = size.text, bty ="n") 

# dev.off()



# ----------------------------------------------
# Abundance indices - predicted vs observed
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

par(mfrow = c(1,1), mar=c(5,5,1,1))

col = c("lightgray","gray")

# x = "I_pred"
x = "I"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
n_plot <- dim(X)[2]

q.mean <- apply(X,2,mean)
q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
  q.inf[k] <- quantile(X[,k],probs=0.05)
  q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
  q.inf1[k] <- quantile(X[,k],probs=0.25)
  q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=Year[1:n_plot],y=q.sup*1.2,ylim=c(0,max(q.sup*1.2)),ylab="Abundance indices",type="n",las=1, main = "Fit Ab Indices")
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup,rev(q.inf)),col = col[1],border=NA)
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=Year[1:n_plot], y=q.mean, col="red", lwd=2, type="l")
points(Year[1:n_plot],data$I_obs[1:n_plot],pch=20,cex=1,col="blue",type="p")

legend(	legend = c("Cred interv. (95%)","IA pred.","IA obs."),
        col = c("lightgray","red","blue"), 	
        lty = c(3,2,NA), pch=c(NA,NA,20), lwd = c(2,2,NA), x = "topright", cex = size.text, bty ="n") 

# dev.off()



#Utilisé le principe du code suivant pour les priors posteriro
# ----------------------------------------------------------------------------
# Marginal posterior distribution - MSY  
# ----------------------------------------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

# res <- 6
# name_figure <- "PostMSY.png"
# png(filename = name_figure, height = 500*res, width = 500*res, res=72*res)

par(mar=c(5,5,1,1), bty="n")

densMSY <- density(mcmc.table$C_MSY)
densMSYp <- density(mcmc.table$C_MSY_p)

plot(densMSY ,ylab = "", xlab = "", xaxt="n", yaxt="n", main="", xlim = c(0,1000), col = "red", type = "l", lty = 1, lwd = 2)
points(densMSYp ,ylab = "", xlab = "", xaxt="n", yaxt="n", main="", xlim = c(0,1000), col = "blue", type = "l", lty = 2, lwd = 2)

axis(side = 1, tick = T, cex.axis=1, las=1)
axis(side = 2, tick = T, lwd.ticks=0, labels=NA, cex.axis=1)
mtext(side=1, "MSY (x 1000 tons)", bty="n", line=3, cex = size.text)
mtext(side=2, "Density", bty="n", line=3, cex = size.text)

legend(	legend = c("Prior","Posterior"),
        col = c("blue","red"), 	
        lty = c(2,1), lwd = c(2,2), x = "topright", cex = size.text, bty ="n") 

par(def.par)

# dev.off()