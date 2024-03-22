# Author: Chris Wymant, Jan-March 2020

library(ggplot2)
library(tidyverse)
library(rriskDistributions)
library(gridExtra)
library(cubature)
library(dplyr)

# Abbreviations:
# serint =  generation time.
# incper = incubation period

# Fast and safe redefinition of integrate.
integrate2 <- function(...) {
  result <- tryCatch({integrate(...)$value},
                     error = function(e) cubintegrate(...)$integral)
  result
}




# General functions
serint <- function(tau, serint.shape, serint.scale) {
  dweibull(x = tau, shape = serint.shape, scale = serint.scale)
}
prob.asymp <- function(tau, incper.shape, incper.scale) {
  1 - pweibull(tau, shape = incper.shape, scale = incper.scale)
} ###Note that it's probability of being pre-symptomatic in symptomatic individuals
prob.symp <- function(tau, incper.shape, incper.scale) {
  pweibull(tau, shape = incper.shape, scale = incper.scale)
}

# Functions specific to this model:


model.gen.beta.s.div.by.R0 <- function(tau, incper.shape, incper.scale,
                                          serint.shape, serint.scale,
                                          P.a, xp) {
  s.of.tau <- prob.symp(tau = tau, incper.shape = incper.shape,
                        incper.scale = incper.scale)  
  serint(tau = tau, serint.shape = serint.shape,
         serint.scale = serint.scale) /
    (xa*P.a +(1 - P.a) * (s.of.tau + xp * (1 - s.of.tau)))
}


###symptomatic contributions
model.gen.beta.sym.tot <- function(tau, incper.shape, incper.scale,
                                   serint.shape, serint.scale,
                                   P.a, xp, R0) {
  (1 - P.a) * R0 * model.gen.beta.s.div.by.R0(tau = tau,
                                                    incper.shape = incper.shape,
                                                    incper.scale = incper.scale,
                                                    serint.shape = serint.shape,
                                                    serint.scale = serint.scale,
                                                    P.a = P.a,
                                                    xp = xp) *
    prob.symp(tau = tau, incper.shape = incper.shape, incper.scale = incper.scale)
}
###Pre-symptomatic contributions
model.gen.beta.presym.tot <- function(tau, incper.shape, incper.scale,
                                      serint.shape, serint.scale,
                                      P.a, xp, R0) {
  xp * (1 - P.a) * R0 * model.gen.beta.s.div.by.R0(tau = tau,
                                                         incper.shape = incper.shape,
                                                         incper.scale = incper.scale,
                                                         serint.shape = serint.shape,
                                                         serint.scale = serint.scale,
                                                         P.a = P.a,
                                                         xp = xp) *
    prob.asymp(tau = tau, incper.shape = incper.shape, incper.scale = incper.scale)
}




model.gen.full.beta.div.by.R0 <-
  function(tau, serint.shape, serint.scale) {
    serint(tau = tau, serint.shape = serint.shape,
           serint.scale = serint.scale) 
  }



model.gen.solve <- function( R0,P.a,  xp, xa,
                            incper.shape, incper.scale, serint.scale,
                            serint.shape, theta.obs,all_constant = all_constant) {
  

  
  integral.of.model.gen.beta.s.div.by.R0 <-
    integrate2(model.gen.beta.s.div.by.R0, lower = 0, upper = Inf,
              incper.scale = incper.scale, incper.shape = incper.shape,
              serint.shape = serint.shape, serint.scale = serint.scale,
              P.a = P.a, xp = xp)
  

  RS <- integrate2(model.gen.beta.sym.tot, lower = 0, upper = Inf,
                  incper.shape = incper.shape,
                  incper.scale = incper.scale,
                  serint.shape = serint.shape,
                  serint.scale = serint.scale,
                  P.a = P.a, xp = xp, R0 = R0)
  RP <- integrate2(model.gen.beta.presym.tot, lower = 0, upper = Inf,
                  incper.shape = incper.shape,
                  incper.scale = incper.scale,
                  serint.shape = serint.shape,
                  serint.scale = serint.scale,
                  P.a = P.a, xp = xp, R0 = R0)
  RA <- R0* P.a * xa * integral.of.model.gen.beta.s.div.by.R0

  
  diff.R0 <- abs(R0 - RA-RS-RP)

  verbose <- FALSE
  if (verbose) {
  
    cat("Diff between R0 and RS + RP +  RA:", diff.R0, "\n")
  
  }
  tolerance <- 1e-3
 
  stopifnot(diff.R0 < tolerance) 
  tau.test1<- seq(from=0, to=13, by=0.05) # x axis 
  df.p.tau.s.loc1.multiply.allcons.div.by.R0<-data.frame(tau = tau.test1, label = "ptau", beta = vapply(
  tau.test1,function(tau) {model.gen.beta.s.div.by.R0(tau = tau,incper.shape = incper.shape, incper.scale = incper.scale,serint.shape = serint.shape, 
    serint.scale = serint.scale,P.a = P.a, xp = xp)}, numeric(1)))

  stopifnot(max(df.p.tau.s.loc1.multiply.allcons.div.by.R0$beta)*R0/all_constant<1) 

  list( R0 = R0,  RA = RA,RS = RS, 
        RP = RP)
  
}
