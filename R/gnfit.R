gnfit<-function(dat, dist, df=NULL, pr=NULL, threshold=NULL){
dat<-as.numeric(dat)
x<-NULL
z<-list()

  op<-par(mfrow = c(1, 2))
if(is.null(pr)){
  if(dist=="gev" | dist=="gpd") stop("Enter Parameters!")
  else if(dist=="t"){
    if(df>2){
    loc<-mean(dat)
    sc<-sqrt((df-2)*var(dat)/df)
    xdat<-(dat-loc)/sc
    prob<-pt(xdat,df)

    qqplot(qt(ppoints(500),df), xdat,main = paste0("Q-Q plot for ", dist,". distribution-","DF:",df),xlab="",ylab="")
    qqline(xdat, distribution = function(p) qt(p,df), probs = c(0.1, 0.6), col ="blue")

    hist(xdat, breaks=20, prob=TRUE, xlab="x-normalized value",main=paste0(dist,". distribution curve over histogram"))
    curve(dt(x,df), col="blue",  add=TRUE, yaxt="n")
    points(xdat, rep(0, length(xdat)))
    }else stop("DF must be > 2")
  }
  else if(dist=="gum"){
   sc<-sqrt(6*var(dat)/pi^2)
   loc<-mean(dat)-0.577*sc
   pr<-c(loc,sc,0)

  prob<-gevf(pr,dat)
  gev.qq(pr,dat)
  gev.his(pr,dat)
      }
else{
    loc<-mean(dat)
    ifelse(dist=="norm",sc<-sd(dat),NA)
    ifelse(dist=="laplace",sc<-sqrt(var(dat)/2),NA)
    ifelse(dist=="logis",sc<-sqrt(3*var(dat)/pi^2),NA)

    prob<-get(paste0("p",dist))(dat,loc,sc)

    qqplot(get(paste0("q",dist))(ppoints(500),loc,sc), dat,main = paste0("Q-Q plot for ", dist,". distribution"),xlab="",ylab="")
    qqline(dat, distribution = function(p) get(paste0("q",dist))(p,loc,sc), probs = c(0.1, 0.6), col ="blue")

    hist(dat, breaks=15, prob=TRUE, xlab="x-normalized value",main=paste0(dist,". distribution curve over histogram"))
    curve(get(paste0("d",dist))(x,loc,sc), col="blue",  add=TRUE, yaxt="n")
    points(dat, rep(0, length(dat)))
  }}else{
    if(dist=="gev"){
      prob<-gevf(pr,dat)
      gev.qq(pr,dat)
      gev.his(pr,dat)
    }else if(dist=="gum"){
      pr[3]<-0
      prob<-gevf(pr,dat)
      gev.qq(pr,dat)
      gev.his(pr,dat)
  }else if(dist=="gpd")
    if(!is.null(threshold)){
      u<-threshold
      dat<-dat[dat>u]
      prob<-gpdf(pr,u,dat)
      gpd.qq(pr,u,dat)
      gpd.his(pr,u,dat)
    }else stop("threshold is missing!")
  }

  n<-length(dat)
  k<-seq(1:n)

  qnor<-qnorm(sort(prob))
  pnor<-pnorm((qnor-mean(qnor))/sd(qnor))
  w<-round((sum((pnor-(2*k-1)/(2*n))^2)+1/(12*n))*(1+0.5/n),4)

  if (w<0.0275) {
    pval<-1-exp(-13.953+775.5*w-12542.61*w^2)
  }
  else if (w<0.051) {
    pval<-1-exp(-5.903+179.546*w-1515.29*w^2)
  }
  else if (w < 0.092) {
    pval <- exp(0.886 - 31.62 * w + 10.897 * w^2)
  }
  else if (w < 1.1) {
    pval <- exp(1.111 - 34.242 * w + 12.832 * w^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
  }
  z$Wpval<-pval

  A<-(-n-sum((2*k-1)*log(pnor)+(2*n+1-2*k)*log(1-pnor))/n)*(1+0.75/n+2.25/n^2)
  A <- round((1 + 0.75/n + 2.25/n^2) * A,4)
  if (A < 0.2) {
    pval<-1-exp(-13.436+101.14*A-223.73*A^2)
  }
  else if (A < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * A - 59.938 * A^2)
  }
  else if (A < 0.6) {
    pval <- exp(0.9177 - 4.279 * A - 1.38 * A^2)
  }
  else if (A < 10) {
    pval <- exp(1.2937 - 5.709 * A + 0.0186 * A^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
  }

  z$Apval<-pval

  z$Cram<-w
  z$Ander<-A
  message("Test of Hypothesis for ", dist, " distribution")
  message("Cramer-von Misses Statistics:  ",z$Cram,"   P-Value:  ", round(z$Wpval,5))
  message("Anderson-Darling Statistics:   ", z$Ander,"   P-Value:  ", round(z$Apval,5))

  class(z)<-"gnfit"
  invisible(z)
}
