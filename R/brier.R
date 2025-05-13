#' Compute Brier Score for Survival Analysis
#'
#' This function calculates
#'
#' @param data.id A dataset for analysis
#' @param survt description
#' @param survp description
#' @param ee description
#' @param s description
#' @param th_v description
#'
#' @return A list containing:
#'   \item{bs}{}
#'
#'
#' @import survival
#' @export

brier3<-function(data.id,survt,survp,s,th_v) {

  nn=nrow(data.id)
  nr=length(survt)
  mt=length(th_v)

  TL1=data.id$TL; TR1=data.id$TR

  delta=ifelse(TR1==Inf,0,1)

  surv0=surv1=rep(0,nn)
  for(i in 1:nn) {
    if(s<survt[1]) surv0[i]=1
    for(j in 2:nr) {
      if(s<survt[j]) {surv0[i]=survp[i,j]; break}
    }
    for(j in 2:nr) {
      if(TL1[i]<survt[j]) {surv1[i]=survp[i,j]; break}
    }
  }

  
  mark2=ff=condG=matrix(0,nn,4)
  for(k in 1:4){
    th=th_v[k]+s
    for(i in 1:nn){
      for(j in 2:nr){
        if(survt[j]<th)   mark2[i,k]=survp[i,j]/surv0[i]
      }
      for(j in 2:nr){
        if(TL1[i]<th&th<TR1[i]&survt[j]<th) {condG[i,k]=survp[i,j]/surv1[i]}
      }
    }
  }

  
  bs<-rep(0,4)
  for(k in 1:4){
    th=th_v[k]+s
    tot=0; tot1=tot2=tot3=tot4=tot5=0; ind1=ind2=ind3=ind4=ind5=0
    for(i in 1:nn){
       if(TR1[i]>s){
      if(s<=TL1[i])                            {
        tot=tot+1
        if(th<=TL1[i])                                 {tot1=tot1+(1-mark2[i,k])^2; ind1=ind1+1}
        else if(TL1[i]<th&th<=TR1[i]&delta[i]==1)      {tot2=tot2+(1-condG[i,k])*(0-mark2[i,k])^2+condG[i,k]*(1-mark2[i,k])^2; ind2=ind2+1}
        else if(TR1[i]<th)                             {tot3=tot3+(0-mark2[i,k])^2 ;ind3=ind3+1}
        else if(TL1[i]<th&delta[i]==0)                 {tot4=tot4+condG[i,k]*(1-mark2[i,k])^2+(1-condG[i,k])*(0-mark2[i,k])^2 ; ind4=ind4+1}
      }
      if(TL1[i]<s&TR1[i]<th&delta[i]==1)               {tot=tot+surv0[i]; tot5=tot5+surv0[i]*(0-mark2[i,k])^2; ind5=ind5+1}
     
    }}
    bs[k]=(tot1+tot2+tot3+tot4+tot5)/tot
    
      }

  return(list("bs"=bs))

}







































