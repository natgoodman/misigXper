#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized data generation functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## compute coverage from sim
cover=function(sim.byn,n,d0,distribution=cq(sim,d2t,meta,bayes,xxx,uuu,norm),mean.prior,sd.prior,
               conf.level) {
  if (missing(distribution)) distribution=if(missing(mean.prior)) 'd2t' else 'meta'
  else distribution=match.arg(distribution,several.ok=TRUE);
  if (any(cq(meta,bayes) %in% distribution)&&(missing(mean.prior)||missing(sd.prior)))
    stop("distribution contains 'meta' or 'bayes' but mean.prior or sd.prior missing");
  cover=do.call(cbind,lapply(n,function(n) {
    sim=sim.byn[[as.character(n)]];
    d.pop=sim$d.pop;
    m=length(d.pop);
    cover=matrix(nrow=length(conf.level),ncol=0);
    if ('sim' %in% distribution) {
      ci=ci_sim(d.pop,d0=d0,conf.level=conf.level);
      cover=cbind(cover,sim=cover1(d.pop,ci,m));
    }
    if ('d2t' %in% distribution) {
      ci=ci_d2t(n=n,d=d0,simplify=F,conf.level=conf.level);
      cover=cbind(cover,d2t=cover1(d.pop,ci,m));
    }
    if ('meta' %in% distribution) {
      meta=meta_normd2t(mean.prior,sd.prior,d0,n);
      ci=ci_norm(mean=meta$mean,sd=meta$sd,d=d0,conf.level=conf.level);
      cover=cbind(cover,meta=cover1(d.pop,ci,m));
    }
    if ('bayes' %in% distribution) {
      ci=ci_bayes(n,d0,mean.prior,sd.prior);
      cover=cbind(cover,bayes=cover1(d.pop,ci,m));
    }
    ## quick hack for d2t posterior probability of d.pop. I think this is what unif is
    if ('xxx' %in% distribution) {
      ci=ci_xxx(n=n,d=d0,simplify=F,conf.level=conf.level);
      cover=cbind(cover,xxx=cover1(d.pop,ci,m));
    }
    ## quick hack for bayes_unif posterior probability of d.pop
    if ('uuu' %in% distribution) {
      ci=ci_uuu(n=n,d0=d0,u=20,simplify=F,conf.level=conf.level);
      cover=cbind(cover,uuu=cover1(d.pop,ci,m));
    }
    ## quick hack for norm posterior probability of d.pop under unif prior
    if ('norm' %in% distribution) {
      ci=ci_norm(mean=d0,sd=sqrt(2/n),conf.level=conf.level);
      cover=cbind(cover,norm=cover1(d.pop,ci,m));
    }
    colnames(cover)=paste(sep='.',colnames(cover),n);
    cover;
  }));
  rownames(cover)=conf.level;
  cover;
}
cover1=function(d.pop,ci,m) {
  lo=ci[1,]; hi=ci[2,];
  cover=do.call(cbind,lapply(d.pop,function(d.pop) between(d.pop,lo,hi)));
  rowSums(cover)/m;
}
## compute concurve from sim. 
concurve=function(sim.byn,n,d0,distribution=cq(sim,d2t,meta,bayes,xxx),mean.prior,sd.prior,
                  conf.level) {
  if (missing(distribution)) distribution=if(missing(mean.prior)) 'd2t' else 'meta'
  else distribution=match.arg(distribution,several.ok=TRUE);
  if (any(cq(meta,bayes) %in% distribution)&&(missing(mean.prior)||missing(sd.prior)))
    stop("distribution contains 'meta' or 'bayes' but mean.prior or sd.prior missing");
  concurve=do.call(cbind,lapply(n,function(n) {
    sim=sim.byn[[as.character(n)]];
    d.pop=sim$d.pop;
    m=nrow(sim);
    concurve=matrix(nrow=length(conf.level),ncol=0);
    if ('sim' %in% distribution) {
      ci=ci_sim(d.pop,d=d0,conf.level=conf.level);
      lo=ci[1,]; hi=ci[2,];
      concurve=cbind(concurve,sim=lo,sim=hi);
    }
    if ('d2t' %in% distribution) {
      ci=ci_d2t(n=n,d=d0,simplify=F,conf.level=conf.level);
      lo=ci[1,]; hi=ci[2,];
      concurve=cbind(concurve,d2t=lo,d2t=hi);
    }
    if ('meta' %in% distribution) {
      meta=meta_normd2t(mean.prior,sd.prior,d0,n);
      ci=ci_norm(mean=meta$mean,sd=meta$sd,d=d0,conf.level=conf.level);
      lo=ci[1,]; hi=ci[2,];
      concurve=cbind(concurve,meta=lo,meta=hi);
    }
    if ('bayes' %in% distribution) {
      ci=ci_bayes(n,d0,mean.prior,sd.prior);
      lo=ci[1,]; hi=ci[2,];
      concurve=cbind(concurve,bayes=lo,bayes=hi);
    }
    ## quick hack for d2t posterior probability of d.pop. I think this is what unif is
    if ('xxx' %in% distribution) {
      ci=ci_xxx(n=n,d=d0,simplify=F,conf.level=conf.level);
      lo=ci[1,]; hi=ci[2,];
      concurve=cbind(concurve,xxx=lo,xxx=hi);
    }
    colnames(concurve)=paste(sep='.',colnames(concurve),n,cq(lo,hi));
    concurve;
  }));
  rownames(concurve)=conf.level;
  concurve;
}
