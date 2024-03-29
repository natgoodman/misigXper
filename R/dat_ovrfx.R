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
## Generate data for ovrfx document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Data Generation for ovrfx ----
dat_ovrfx=function() {
  ## do simulation
  param(n.rand,m.rand,d.rand);
  dosim_rand(n.rand,m.rand,d.rand);
  param(n.fixd,m.fixd,d.fixd);
  dosim_fixd(n.fixd,m.fixd,d.fixd);
  ## compute mean significant effect size
  domeand_fixd(n.fixd,d.fixd);
  domeand_d2t(n.fixd,d.fixd);
  invisible();
}
