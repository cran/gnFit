rskFac<-function(dat,alpha=0.1,dist="norm",df=NULL)
{
  z<-list()
  qnt<-function(alpha=0.1,dist="norm",df=NULL){
    z<-list()

    if(dist=="t"){
      z$c1<-1
      z$c2<-sqrt((df-2)/df)
      z$qn<-qt(alpha,df)
    }
    else if (dist=="norm"){
      z$c1<-1
      z$c2<-1
      z$qn<-qnorm(alpha)
    }
    else if (dist=="logis"){
      z$c1<-1
      z$c2<-sqrt(3)/pi
      z$qn<-qlogis(alpha)
    }
    else if (dist=="laplace"){
      z$c1<-1
      z$c2<-sqrt(1/2)
      z$qn<-qlaplace(alpha)
    }
    else if (dist=="gum"){
      #z$c1<-1-0.577*pi^2(sqrt(6)*retn*std)
      z$c1<-1
      z$c2<-sqrt(6)/pi
      z$qn<-gum.q(1-alpha,0,1)
    }

    class(z)<-"qnt"
    invisible(z)
  }

  dat.port<-as.numeric(dat)
  ret<-mean(dat.port)
  std<-sqrt(var(dat.port))
  if(dist=="t" & is.null(df)) stop("DF is missing!")

  #dat<-(na.omit(RRSC(p,sub,L,D)$return))
  #z$dat.port<-t((w0)%*%t(dat))

  #  dat.port<-datPort(p,w0,sub,L,D)$datport

  #  if(length(p)==1) w0<-1
  #  else w0<-t(w0)/100



 # try(gnfit(as.numeric(dat.port),dist,df),silent=TRUE)


  qnt0<-qnt(alpha,dist,df)
  qnt1<-qnt(1-alpha,dist,df)

  if(dist=="gum")
  {
    expVal<-qnt0$c1*ret-0.577*sqrt(6)*std/(pi)
    #stdVal<-std*qnt0$c2
    VaR<-round((qnt0$c1*ret-0.577*sqrt(6)*std/(pi)+qnt0$qn*std*qnt0$c2),3)
    VaR1<-round((qnt1$c1*ret-0.577*sqrt(6)*std/(pi)-qnt1$qn*std*qnt1$c2),3)
    #VaRP<-round(-qnt0$c1*ret-0.577*sqrt(6)*std/(pi)-qnt0$qn*std*qnt0$c2,3)
  }
  else
  {

    #  qnt0<-qnt(alpha,dist,df)
    expVal<-qnt0$c1*ret
    #stdVal<-std*qnt0$c2
    VaR<-round((qnt0$c1*ret+qnt0$qn*std*qnt0$c2),3)
    VaR1<-round((qnt1$c1*ret-qnt1$qn*std*qnt1$c2),3)
#    VaR1<-round((qnt1$c1*ret+qnt1$qn*std*qnt1$c2),3)
    #  qnt01<-qnt(1-alpha,dist,df)
    #  VaRP<-round((qnt01$c1*ret+qnt01$qn*std*qnt01$c2),3)

  }

  VaR<-as.numeric(VaR)
  VaR1<-as.numeric(VaR1)
  #VaR<--0.08
  #dat.port<-na.omit(RRSC("FRVR1",sub="2016::")$return)

d<-pmax(VaR-dat.port,0)
AVaR_N<-round(VaR-mean(d)/(alpha),4)


  #dP<-pmax(VaR+dat.port,0)
  #AVaR_P<-round(-VaR+mean(dP)/(1-alpha),4)

d<-pmax(VaR+dat.port,0)
AVaR_P<-round(-VaR+mean(d)/(1-alpha),4)

  #AVaR_P3<-round((expVal-(alpha)*AVaR_N)/(1-alpha),4)
  #AVaR_P4<-round((expVal-(alpha)*AVaR_P)/(1-alpha),4)
  #E<-round(-(1-alpha)*AVaR_P+alpha*AVaR_N,4)
  #VaRD<-round(expVal-VaR,4)
  message("Under assumption *", dist,"* df,")
  message("E(X) = ",round(expVal,4)*100," (%) ")
  message("VaR(X) = ", VaR*100," (%), at level alpha = ", alpha )

  #message("VaRD(X) = ", VaRD," at level (1-alpha) = ", 1-alpha )
  #message("  VaR at lower tail: ",VaR)
  message("Lower AVaR(X) = ",round(AVaR_N*100,2),"(%), at level alpha = ", alpha)
  #message("Upper tail at level alpha= ", alpha )
#  message("  VaR at upper tail: ",-VaR1)
  message("Upper AVaR(X) =  ",round(AVaR_P*100,2),"(%), at level alpha = ", alpha)

  #message(" AVaR(X) = ",round(AVaR_P3*100,2),"     (%), at level 1-alpha = ", 1-alpha)
  #message(" AVaR(X) = ",round(AVaR_P4*100,2),"     (%), at level 1-alpha = ", 1-alpha)
# message("AVaRD(X) = ",AVaR_P2," at level (1-alpha) = ", 1-alpha)
  #message("Expected VaR: ",E,"\n" )

  #message("SR of Return Port. (%): ",round(ret/std*100,4))
  #message("Volatility of Port. (%): ",round(std*100,4),"\n" )

  #  }
  z$VaR<-VaR
  z$AVaR_p<-AVaR_P
  z$AVaR_n<-AVaR_N
  #z$EVaR<-E
  z$mRtn<-round(ret,4)
  class(z)<-"rskFac"
  invisible(z)
}
