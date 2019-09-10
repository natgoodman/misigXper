#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
library(nor1mix);
## TODO: turn constants into params
dat_confi=
  function(n=param(n.confi),m=3000,d0=param(d0.confi),
           mean.true=0.3,mean.false=0,sd.true=0.1,sd.false=0.05,
           n.lo=c(10,20),n.hi=200,prop.lo=0.25,prop.hi=0.75,...) {
    ## do normals
    dosim(n=n.lo[1],m=m,d0=d0,id='lo_norm',d.gen=rnorm,
          d.args=list(mean==mean.true,sd=sd.true),...);
    dosim(n=n.lo[2],m=m,d0=d0,id='lo_norm',d.gen=rnorm,
          d.args=list(mean==mean.true,sd=sd.true),...);
    dosim(n=n.hi,m=m,d0=d0,id='hi_norm',d.gen=rnorm,
          d.args=list(mean==mean.true,sd=sd.true),...);
    ## do mixtures
    mean.mix=c(mean.true,mean.false);
    sd.mix=c(sd.true,sd.false);
    dosim(n=n.lo[1],m=m,d0=d0,id='lo_lo',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.lo,1-prop.lo))),...);
    dosim(n=n.lo[2],m=m,d0=d0,id='lo_lo',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.lo,1-prop.lo))),...);
    dosim(n=n.hi,m=m,d0=d0,id='hi_lo',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.lo,1-prop.lo))),...);
    dosim(n=n.lo[1],m=m,d0=d0,id='lo_hi',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.hi,1-prop.hi))),...);
    dosim(n=n.lo[2],m=m,d0=d0,id='lo_hi',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.hi,1-prop.hi))),...);
    dosim(n=n.hi,m=m,d0=d0,id='hi_hi',d.gen=rnorMix,
          d.args=list(obj=norMix(mu=mean.mix,sigma=sd.mix,w=c(prop.hi,1-prop.hi))),...);
    invisible();
  }
vrnorm=Vectorize(rnorm,"mean");

## NG 19-09-04: finally broke down and rewrote dosim to do pruning and hit target number of rows
## TODO: turn constants into params
dosim=
  function(n,m=3000,d0=0.5,id='mix',d.gen=rnorMix,d.args=list(),tol=1e-3,m1=1e4,mmax=1e8) {
    param(verbose,debug);
    if (verbose) {
      if (identical(d.gen,rnorMix))
        expect=paste('expect',
                     round(emix(mix=d.args$obj,n=n,m=m,d0=d0,tol=tol,mmax=mmax)/m1),'iters')
      else expect=NULL;
      print(paste(
        collapse=' ',c('>>> dosim:',nvq(id,n,m,d0),expect,format(Sys.time(),"%b %d %X"))));
    }
    sim=data.frame(row.names=NULL);
    m1.sum=0; i=0;
    while(nrow(sim)<m&&m1.sum<mmax) {
      m1=min(m1,mmax-m1.sum);           # overly cautious, but why not?
      group0=replicate(m1,rnorm(n,mean=0));
      d=do.call(d.gen,c(n=m1,d.args));
      group1=vrnorm(n,mean=d);
      mean0=colMeans(group0);
      mean1=colMeans(group1);
      d.raw=mean1-mean0;
      sd0=apply(group0,2,sd);
      sd1=apply(group1,2,sd);
      sd=pooled_sd(sd0,sd1);
      d.sdz=d.raw/sd;
      sim1=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,row.names=NULL);
      sim=rbind(sim,subset(sim1,subset=near(d.sdz,d0,tol)));
      m1.sum=m1.sum+m1; 
      if (debug) {
        i=i+1;
        print(paste(sep=' ','+++ dosim:',nvq(i),paste_nv('nrow',nrow(sim)),expect));
      }
    }
    if (nrow(sim)<m)
      warning(paste('dosim failed to generate enough rows. wanted',m,'got',nrow(sim)))
    else sim=sim[1:m,];
    save_sim(sim,n,m,d0,id);
    invisible(sim);
  }

## file functions. should be in datman but I'm worried about name conflicts...
filename_sim=function(n,m,d0,id) 
  filename(param(simdir),base=paste(sep='_','sim',id),
           tail=paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
           suffix='RData');
save_sim=function(sim,n,m,d0,id) save(sim,file=filename_sim(n,m,d0,id));

load_sim=get_sim=function(n,m,d0,id,tol=param(tol),prune=F) {
  sim=load_(filename_sim(n,m,d0,id),'sim')
  if (prune) sim=subset(sim,subset=near(d.sdz,d0,tol));
  invisible(sim);
}

## for debug output
emix=function(mix,n,m,d0,tol=1e-3,mmax=1e8) {
  f=function()
    integrate(function(d.pop) dnorMix(d.pop,mix)*d_d2t(n=n,d=d0,d0=d.pop),-Inf,Inf)$value*2*tol;
  uniroot(function(m.need) m.need*f()-m,interval=c(1,mmax))$root;
}
