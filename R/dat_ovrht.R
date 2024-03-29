#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-19
##          from dodata.R created 19-02-18
##          from ovrfx.R created 19-02-03 
##          from siglo.R 19-01-01
##          from repwr/R/repwr.R created 17-10-05 
##           and repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for ovrht document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Data Generation for ovrht ----
dat_ovrht=function() {
  param(n.hetd,m.hetd,d.hetd,sd.hetd);
  dosim_hetd(n.hetd,m.hetd,d.hetd,sd.hetd);
  ## generate pval and ci tables
  param(n.ovrht,d.ovrht,sd.ovrht,sig.dat);
  dopval_d2ht(n.ovrht,sd.ovrht,sig.dat);
  doci_d2ht(n.ovrht,d.ovrht,sd.ovrht);
  invisible();
}
