#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Sandbox program for confidence intervals.
## NOT USED in running code but has code I don't want to lose just yet
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
source('R/source.R');                   # source the others

## ---- confi ----

## ---- Simulation Functions ----
dat_confi=
  ## function(n=c(10,seq(20,200,by=20),400),m=1e6,d0=0.5,
  ##          id='unif',d.gen=runif,d.args=list(min=d0-3,max=d0+3),
  ##          verbose=TRUE) {
  ## function(n=c(10,seq(20,200,by=20),400),m=1e6,d0=0.5,verbose=TRUE,...) {
  function(n=c(20,60,100,200),m=1e6,d0=0.5,verbose=TRUE,...) {
    param(n=n,m=m,d0=d0,verbose=verbose,doc='confi',datadir='data/confi');
    cases=expand.grid(n=n,m=m,d0=d0);
    withrows(cases,case,dosim(n=n,m=m,d0=d0,verbose=verbose,...));
}
vrnorm=Vectorize(rnorm,"mean");
dosim=
  function(n,m,d0=0.5,id='unif',d.gen=runif,d.args=list(min=d0-3,max=d0+3),m0=1e4,verbose=T) {
    more=m; i=1;
    while(more>0) {
      m1=min(m0,more);
      if (verbose) print(paste(sep=' ','>>> dosim:',nvq(id,i,n,m,d0,more,m1)));
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
      sim=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,row.names=NULL);
      save_sim(sim,n,m,d0,id,i);
      more=more-m1; i=i+1;
    }
    ## consolidate subfiles into one
    sim=cat_sim(n,m,d0,id);
    invisible(sim);
  }
doc_confi=
  function(n=param(n),m=param(m),d0=param(d0),id='unif',tol=1e-3,smooth='aspline',
           conf.level=seq(0.05,0.95,by=0.05)) {
    init_doc(figscreen=T);
    cases=expand.grid(n=n,m=m,d0=d0);
    sim=do.call(rbind,withrows(cases,case,{
      sim=load_sim(n=n,m=m,d0=d0,id=id);
      sim=subset(sim,subset=near(d.sdz,d0,tol));
    }));
    sim.byn=split(sim,sim$n);
    sim<<-sim;                          # for debugging
    sim.byn<<-sim.byn;                  # for debugging
    cover=do.call(cbind,lapply(n,function(n) {
      sim=sim.byn[[as.character(n)]];
      sapply(conf.level,function(conf.level) {
        ci=ci_d2t(n=n,d=sim$d.sdz,simplify=F,conf.level=conf.level);
        cover=between(sim$d.pop,ci[1,],ci[2,]);
        sum(cover)/nrow(sim);
      })}));
    colnames(cover)=names(sim.byn);
    rownames(cover)=conf.level;
    cover<<-cover;                      # for debugging
    title=figtitle('coverage vs. confidence level by n',id=id,d0=d_pretty(d0),m=m_pretty(m));
    figname=paste(sep='.',
                  'cover_by_n',id,
                  paste(sep=',',paste_nv(d0,d_pretty(d0)),paste_nv(m,m_pretty(m))));
    dofig(plotm_cover,figname,cover=cover,title=title);
    cover;
  }
## figure 1 (some version...). histogram, d2t sampling distribution, ci's
## legend tells whether to draw pval legend
##   legend.x0 is harcoded x-position of legend - use with care!
##   legend.xscale, legend.yscale are fractions of plat areas
##   legend.cex is cex for legend text
## col.d2t, etc are vectors: 1st elements for sig, 2nd for nonsig
##   if NULL, col.d2t means color by pval
## col.ci, lty.ci, lwd.ci are vectors: 1st elements for ci's that cover d0, 2nd ones that don't
## Note: xlim, ylim often passed in ...
plothist_ci=plothist_confx=
  function(sim,n=sim$n[1],d.pop=sim$d.pop[1],conf.level=0.95,
           p=seq(0.025,0.975,by=0.025),d=q_d2t(n=n,d0=d.pop,p=p),
           title,xlab="observed effect size (Cohen's d)",ylab='probability density',
           legend=T,legend.xscale=1/8,legend.yscale=1/3,legend.x0=NULL,legend.cex=0.75,
           hist.col='grey90',border.col='grey80',
           d2t.col=NULL,d2t.lty='solid',d2t.lwd=c(4,2),
           ci.col=cq(green,grey50),ci.lty=cq(solid,dotted),ci.lwd=c(0.5,2),
           ci.pch=19,ci.cex=c(0.5,0.75),ci.tol=1e-5,
           vline=cq(d.pop,d.conf),
           ...) {
    d2t.lty=rep(d2t.lty,len=2); d2t.lwd=rep(d2t.lwd,len=2);
    ci.col=rep(ci.col,len=2); ci.lty=rep(ci.lty,len=2); ci.lwd=rep(ci.lwd,len=2);
    ci.pch=rep(ci.pch,len=2); ci.cex=rep(ci.cex,len=2);
    
    plothist(sim=sim,col=hist.col,border=border.col,title=title,xlab=xlab,ylab=ylab,legend=F,...);
    plotpvsd(add=T,n=n,d0=d.pop,col=d2t.col,
             lty.sig=d2t.lty[1],lty.nonsig=d2t.lty[2],lwd.sig=d2t.lwd[1],lwd.nonsig=d2t.lwd[2]);
    if (legend)
      pval_legend(x.scale=legend.xscale,y.scale=legend.yscale,x0=legend.x0,cex=legend.cex);
    d.conf=d_conf(n=n,d0=d.pop,conf.level=conf.level);
    ## if (missing(d)) {
    ##   ## TODO: extend wings until density too low
    ##   d.diff=diff(d.conf)/10;
    ##   d=c(d.conf[1]-d.diff,seq(d.conf[1],d.conf[2],len=10),d.conf[2]+d.diff);
    ## }
    ci=ci_d2t(n=n,d=d,conf.level=conf.level);
    y=d_d2t(n=n,d=d,d0=d.pop);
    col=ifelse(between(d.pop,ci[1,],ci[2,],tol),ci.col[1],ci.col[2]);
    lty=ifelse(between(d.pop,ci[1,],ci[2,],tol),ci.lty[1],ci.lty[2]);
    lwd=ifelse(between(d.pop,ci[1,],ci[2,],tol),ci.lwd[1],ci.lwd[2]);
    pch=ifelse(between(d.pop,ci[1,],ci[2,],tol),ci.pch[1],ci.pch[2]);
    cex=ifelse(between(d.pop,ci[1,],ci[2,],tol),ci.cex[1],ci.cex[2]);
    hline(y=y,x0=ci[1,],x=ci[2,],col=col,lty=lty,lwd=lwd);
    points(d,y,col=col,pch=pch,cex=cex);
    if ('d.pop' %in% vline) vhline(vline=d.pop,col='black',lty='solid');
    if ('d.conf' %in% vline) vhline(vline=d.conf,col=ci.col[1],lty='dotted');
  }
d_conf=function(n,d0,simplify=T,conf.level=param(conf.level)) {
  d.conf=d_conf_(n,d0,conf.level);
  if (simplify&ncol(d.conf)==1) d.conf=as.vector(d.conf);
  d.conf;
}
d_conf_=Vectorize(function(n,d0,conf.level) {
  lo=uniroot(function(d) ci_d2t(n=n,d=d,conf.level=conf.level)[2]-d0,interval=c(-5,d0))$root;
  hi=uniroot(function(d) ci_d2t(n=n,d=d,conf.level=conf.level)[1]-d0,interval=c(d0,5))$root;
  c(lo,hi);
})
## figure 1 (some version...).
## histogram, d2t sampling distribution (in grey), quantiles
## Note: xlim, ylim often passed in ...
plothist_q=
  function(sim,n=sim$n[1],d.pop=sim$d.pop[1],quantile=c(0.5,0.8,0.95),
           title,xlab="observed effect size (Cohen's d)",ylab='probability density',
           ## legend=T,legend.xscale=1/8,legend.yscale=1/3,legend.x0=NULL,legend.cex=0.75,
           hist.col='grey90',border.col='grey80',
           distr.col='black',distr.lty='solid',distr.lwd=2,
           q.col='black',q.lty='dotted',q.lwd=2,
           dpop.col='black',dpop.lty='solid',dpop.lwd=1,
           ...) {
    plothist(sim=sim,col=hist.col,border=border.col,title=title,xlab=xlab,ylab=ylab,legend=F,...);
    f=function(x) d_d2t(n,d=x,d0=d.pop);
    curve(f,add=T,col=distr.col,lty=distr.lty,lwd=distr.lwd);
    ## draw quantiles
    p=(1-quantile)/2;
    x=q_d2t(n=n,d0=d.pop,p=p,lower.tail=TRUE); y=d_d2t(n=n,d0=d.pop,d=x);
    vline(x=x,y=y,col=q.col,lty=q.lty,lwd=q.lwd);
    label=paste0(100*quantile,'% ');
    text(x,y,label,col=q.col,lty=q.lty,lwd=q.lwd,cex=0.75,adj=c(1,0));
    x=q_d2t(n=n,d0=d.pop,p=p,lower.tail=FALSE); y=d_d2t(n=n,d0=d.pop,d=x);
    vline(x=x,y=y,col=q.col,lty=q.lty,lwd=q.lwd);
    label=paste0(' ',100*quantile,'%');
    text(x,y,label,col=q.col,lty=q.lty,lwd=q.lwd,cex=0.75,adj=c(0,0));
    ## draw d.pop
    ## vhline(vline=d.pop,col=dpop.col,lty=dpop.lty,lwd=dpop.lwd);
  }



## wrapper for plotm to show coverage from simulated data
plotm_cover=
  function(cover,title=NULL,smooth='aspline',
           col=RColorBrewer::brewer.pal(max(3,ncol(cover)),'Set1'),lty='solid',lwd=2,
           legend='topleft',xlab='confidence level',ylab='coverage',
           col.diag='grey',lty.diag='dotted',lwd.diag=1.5) {
    x=as.numeric(rownames(cover));
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
    legend.labels=paste(sep='=','n',colnames(cover));
    plotm(x=x,y=cover,title=title,col=col,lty=lty,lwd=lwd,smooth=smooth,
          legend=legend,legend.labels=legend.labels,xlab=xlab,ylab=ylab);
    abline(a=0,b=1,col=col.diag,lty=lty.diag,lwd=lwd.diag);
    grid();
  }

filename_sim=function(n,m,d0,id,i=NULL) 
  filename(param(datadir),'sim',base=paste(sep='_','sim',id),
           tail=paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
           suffix=paste(collapse='.',c(if (!is.null(i)) sprintf("%03i",i),'RData')));

save_sim=function(sim,n,m,d0,id,i=NULL) {
  save(sim,file=filename_sim(n,m,d0,id,i));
}
load_sim=get_sim=function(n,m,d0,id,i=NULL,subfiles=FALSE) {
  if ((is.null(i)&!subfiles)|!is.null(i)) sim=load_(filename_sim(n,m,d0,id,i),'sim')
  ## if (!is.null(i)) sim=load_(filename_sim(n,m,d0,id,i),'sim')
  else {
     pattern=paste0('sim_',id,'\\.',
                    paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
                    paste0('\\.','\\d+','\\.','\\RData'));
    files=
      list.files('data/confi/sim',full.names=T,pattern=pattern);
    if (length(files)==0) stop(paste0('no files found: pattern=',pattern));
    sim=do.call(rbind,lapply(files,function(file) load_(file,'sim')))
  }
  invisible(sim);
}
## remove sim subfiles after consolidating
rm_sim=function(n,m,d0,id,i=NULL,subfiles=T) {
  if ((is.null(i)&!subfiles)|!is.null(i)) file.remove(filename_sim(n,m,d0,id,i))
       pattern=paste0('sim_',id,'\\.',
                    paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
                    paste0('\\.','\\d+','\\.','\\RData'));
  files=
      list.files('data/confi/sim',full.names=T,pattern=pattern);
  if (length(files)==0) stop(paste0('no files found: pattern=',pattern));
  file.remove(files);
}
## consolidate sim subfiles, save as one file, rm subfiles
cat_sim=function(n,m,d0,id) {
  sim=load_sim(n,m,d0,id,subfiles=T);
  save_sim(sim,n,m,d0,id);
  rm_sim(n,m,d0,id);
  invisible(sim);
}
## load all (consolidated) sim files for parameter grid, optionally pruning to region near d0
load_sim_all=get_sim_all=
  function(n=param(n.confi),m=param(m.confi),d0=param(d0.confi),id=param(id.confi),
           tol=param(tol),prune=T) {
    cases=expand.grid(n=n,m=m,d0=d0);
    sim=do.call(rbind,withrows(cases,case,{
      sim=load_sim(n=n,m=m,d0=d0,id=id);
      if (prune) sim=subset(sim,subset=near(d.sdz,d0,tol));
    }));
    invisible(sim);
  }
