smry.t.test=function(n,mean0,mean1,sd0,sd1,conf.level=0.95) {
  mx=mean0; my=mean1;
  df <- 2*(n-1)
  stderr=sqrt((sd0^2+sd1^2)/n)
  tstat <- (mx - my)/stderr
  pval <- 2 * pt(-abs(tstat), df)
  alpha <- 1 - conf.level
  cint <- qt(1 - alpha/2, df)
  cint <- tstat + c(-cint, cint)
  
  cint <- cint * stderr
  rval <- list(statistic = tstat, parameter = df, p.value = pval, conf.int = cint)
  rval
}

## ci_tt=function(n,mean0,mean1,sd0,sd1,conf.level=param(conf.level)) {
##   df=2*(n-1)
##   stderr=sqrt((sd0^2+sd1^2)/n)
##   tstat <- (mean0 - mean1)/stderr
##   pval <- 2 * pt(-abs(tstat), df)
##   alpha <- 1 - conf.level
##   cint <- qt(1 - alpha/2, df)
##   cint <- tstat + c(-cint, cint)
##   cint <- cint * stderr
##   cint;
## }

ci_tt=function(n,d.raw,sd0,sd1,simplify=T,conf.level=param(conf.level)) {
  ci=ci_tt_(n,d.raw,sd0,sd1,conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}

ci_tt_=Vectorize(function(n,d.raw,sd0,sd1,conf.level) {
  df=2*(n-1);
  stderr=sqrt((sd0^2+sd1^2)/n);
  tstat=d.raw/stderr;
  alpha=1-conf.level;
  ci=qt(1-alpha/2,df);
  lo=(tstat-ci)*stderr;
  hi=(tstat+ci)*stderr;
  c(lo,hi);
})


posterior=function(d.true,n,d.obs,prior) {
  ## prior probability of d.true
  pA=prior;
  ## likelihood - probability of d.obs given d.true
  pBgivenA=function(d.true) d_d2t(n=n,d0=d.true,d=d.obs);
  ## joint prob density
  pAgivenB=function(d.true) pA(d.true)*pBgivenA(d.true);
  ## average liklihood - denominator in Bayes formula
  pB=integrate(function(d.true) pAgivenB(d.true),-Inf,Inf)$value;
  ## final answer
  pAgivenB(d.true)/pB;
}
posterior=function(d.true,n,d.obs,prior) {
  ## probability of d.obs given d.true for the case at hand:
  ##   two group difference of means, equal size and variance
  P_obsGIVENtrue=function(d.true) d_d2t(n=n,d0=d.true,d=d.obs);
  ## numerator in Bayes formula
  numerator=function(d.true) prior(d.true)*P_obsGIVENtrue(d.true);
  ## denominator in Bayes formula
  P_obs=integrate(function(d.true) numerator(d.true),-Inf,Inf)$value;
  ## final answer
  numerator(d.true)/P_obs;
}

