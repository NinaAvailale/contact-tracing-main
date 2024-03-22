#####################################
# This code solves for the different 
# functional forms contributing to from specified
# beta(tau) in our model of infectiousness.

# (Ferretti, Wymant et al, Science 2020)
#####################################


###Based on this,We further calculated the transmission dynamics in various contact locations###



library(ggplot2)
library(tidyverse)
library(rriskDistributions)
library(gridExtra)
library(RColorBrewer)
library(R.utils)


# Abbreviations:
# serint = serial interval. Deprecated. It now in fact refers to generation time.
# incper = incubation period

# TODO: set this to the directory you want to work in. Make sure
# Theta_EstimatedByRenewalEquation_funcs.R is in there.
ncov.dir <- "C:\\contact_tracing_main\\model" 

setwd(ncov.dir)
source("Theta_EstimatedByRenewalEquation_funcs_my_fig2.R")

# A tibble for parameters - their central values, and lower and upper CIs (for
# when we assign CIs to our calculations.)
data.params <- tibble(name = character(0),
                      central = numeric(0),
                      lower   = numeric(0),
                      upper   = numeric(0))
################################################################################
# INPUT DATA
# Add each parameter and 95% CIs as a new row to the data.params data frame.



# Add Incubation period parameters for weibull distribution
# All have units of days
data.params <- bind_rows(data.params,
                         list(name = "incper.shape", central = 1.792,
                              lower = 1.719, upper = 1.856)
)
data.params <- bind_rows(data.params,
                         list(name = "incper.scale", central = 3.881,
                              lower = 3.817, upper = 3.984)
)
# Add Generation time parameters for weibull distribution
# All have units of days
data.params <- bind_rows(data.params,
                         list(name = "serint.shape", central = 1.364, # Weibull 
                              lower = 1.263, upper = 1.469)
) 
data.params <- bind_rows(data.params,
                         list(name = "serint.scale",  central = 3.233, 
                              lower = 2.993, upper = 3.471) 
) 

# Input Known Params  
xp <- 1;
P.a <- 0.218; 
xa <- 1;R0<-10

# Plotting params:
tau.test <- seq(from=0, to=13, by=0.05) # x axis 
colors <- brewer.pal(4, "Paired")
names(colors) <- c("pre-symptomatic", "symptomatic", "environmental", "asymptomatic")
base.font.size <- 27
plot.height <- 8
plot.width <- 8
minor_breaks <- seq(0, 13, 2)
breaks <- seq(0, 13, 2)

################################################################################

# Get the median of the serint distribution. 
serint.median <- qweibull(p = 0.5, scale =
                          data.params[data.params$name == "serint.scale", ]$central,
                          shape = data.params[data.params$name == "serint.shape", ]$central)

plot(tau.test, serint(tau = tau.test,
                     serint.shape = data.params[data.params$name == "serint.shape", ]$central,
                     serint.scale = data.params[data.params$name == "serint.scale", ]$central))


  
  # For convenience, define variables e.g. called foo instead of
  # pararams[data.params$name == "foo", ]$central
  for (param.num in seq(from=1, to=nrow(data.params))) {
    assign(data.params$name[[param.num]], data.params$central[[param.num]])
  }
  

  #1-5:Dwelling,Work,Cospace-time,Community,Other(Unknown)
  ####average daily contact numbers for each contact pattern (derived from contact tracing data collected)
  c1<-2.95;c2<-13.92;c3<-29.08;c4<-25.92; c5<-3.62
  ####average transmission rate constant vs home for each contact patterns (derived from epidemiological data and contact tracing collected)
  x11<-1;x21<-0.0958;x31<-0.0283;x41<-0.0419;x51<-0.3399
  
  
  #####fraction of R0 for each contact pattern
  all_constant<-c1*x11+c2*x21+c3*x31+c4*x41+c5*x51
  frac_c1<-c1*x11/all_constant
  frac_c2<-c2*x21/all_constant
  frac_c3<-c3*x31/all_constant
  frac_c4<-c4*x41/all_constant
  frac_c5<-c5*x51/all_constant

  # Solve!
  dummy <- model.gen.solve(R0 = R0,P.a = P.a,
                           xp = xp, xa = xa, incper.shape = incper.shape,
                           incper.scale = incper.scale,
                           serint.scale = serint.scale,
                           serint.shape = serint.shape,
                           all_constant = all_constant)
  R0 <- dummy$R0
  RA <- dummy$RA # contribution to R0 from asymps
  RS <- dummy$RS # contribution to R0 from symps
  RP <- dummy$RP # contribution to R0 from pre-symps

  
  # Plot the different contributions to beta(tau):
  
  df.beta.p <- data.frame(tau = tau.test, label = "pre-symptomatic", beta = vapply(
    tau.test, model.gen.beta.presym.tot, numeric(1), incper.shape = incper.shape,
    incper.scale = incper.scale, serint.shape = serint.shape,
    serint.scale = serint.scale, P.a = P.a, xp = xp, R0 = R0))
  
  df.beta.s <- data.frame(tau = tau.test, label = "symptomatic", beta = vapply(
    tau.test, model.gen.beta.sym.tot, numeric(1), incper.shape = incper.shape,
    incper.scale = incper.scale, serint.shape = serint.shape,
    serint.scale = serint.scale, P.a = P.a, xp = xp, R0 = R0))
  
  df.beta.a <- data.frame(tau = tau.test, label = "asymptomatic", beta = vapply(
    tau.test, function(tau) {R0 * P.a * xa *model.gen.beta.s.div.by.R0(tau = tau,
    incper.shape = incper.shape,incper.scale = incper.scale,serint.shape = serint.shape,
    serint.scale = serint.scale,P.a = P.a, xp = xp)} , numeric(1)))
  
  df.beta.p[which.max(df.beta.p$beta),1]

  df.plot <- rbind(df.beta.a, df.beta.p, df.beta.s)
  df.plot$label <- factor(df.plot$label, levels = c("pre-symptomatic",
                                                    "symptomatic",
                                                    "asymptomatic"))
  
  # Go from long to wide format for ease of finding the maximum stacked value
  df.plot.wide <- spread(df.plot, label, beta)
  df.plot.wide$max <- df.plot.wide$`pre-symptomatic` + df.plot.wide$symptomatic
                    + df.plot.wide$asymptomatic
  ymax <- max(df.plot.wide$max * 1.5)

  p <- ggplot(df.plot, aes(x=tau, y=beta)) + 
    theme_bw(base_size = base.font.size) +
    geom_area(aes(fill=label)) +
    labs(x = expression(paste(tau, " (days)")),
         y = expression(paste(beta, "(", tau, ")  (transmissions per day)")),
         fill = bquote(paste('R'['0'] * ' = ' * .(format(round(R0, 1), nsmall = 1)) * ":" )),
         title = NULL,  
         subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
    scale_fill_manual(values = colors[c("pre-symptomatic", "symptomatic",  "asymptomatic")],
                      labels = c(
                        bquote(paste('R'['p'] * ' = ' * .(format(round(RP, 1), nsmall = 1)) * " from pre-symptomatic")),
                        bquote(paste('R'['s'] * ' = ' * .(format(round(RS, 1), nsmall = 1)) * " from symptomatic")),
                        bquote(paste('R'['a'] * ' = ' * .(format(round(RA, 1), nsmall = 1)) * " from asymptomatic")))) +
    scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks)
  p
  plot.name <- paste0("supple_figure_for_different_symptom_type.pdf")
  ggsave(plot.name, p, height = plot.height, width = 13)
  # 

  #####transmission for each contact location
  R_c1<-frac_c1*R0
  R_c2<-frac_c2*R0
  R_c3<-frac_c3*R0
  R_c4<-frac_c4*R0
  R_c5<-frac_c5*R0


#######################################################################
#########Plot the different contributions for each contact location to beta(tau):
  
df.p.tau.multiply.allcons.div.by.R0<-data.frame(tau = tau.test, label = "ptau", beta = vapply(
    tau.test,function(tau) {model.gen.full.beta.div.by.R0(tau = tau, serint.shape = serint.shape, serint.scale = serint.scale)}, numeric(1)))
 

  df.p.tau.loc.c1<-data.frame(tau = tau.test,label = "home",beta = df.p.tau.multiply.allcons.div.by.R0$beta*R0/all_constant)  
  df.p.tau.loc.c2<-data.frame(tau = tau.test,label = "work",beta = df.p.tau.loc.c1$beta*x21)
  df.p.tau.loc.c3<-data.frame(tau = tau.test,label = "social",beta = df.p.tau.loc.c1$beta*x31)
  df.p.tau.loc.c4<-data.frame(tau = tau.test,label = "community",beta = df.p.tau.loc.c1$beta*x41)
  df.p.tau.loc.c5<-data.frame(tau = tau.test,label = "other",beta = df.p.tau.loc.c1$beta*x51)


  df.plot2 <- rbind(df.p.tau.loc.c1, df.p.tau.loc.c2, df.p.tau.loc.c3,df.p.tau.loc.c4, df.p.tau.loc.c5)
  df.plot2$label <- factor(df.plot2$label, levels = c("home",
                                                      "work",
                                                      "social",
                                                      "community",
                                                      "other"))
# Go from long to wide format for ease of finding the maximum stacked value
df.plot.wide2 <- spread(df.plot2, label, beta)
df.plot.wide2$max <- df.plot.wide2$home + df.plot.wide2$work+ df.plot.wide2$social
+ df.plot.wide2$community + df.plot.wide2$other
ymax2 <- max(df.plot.wide2$max * 1.5)
colors <- brewer.pal(8, "Set2")
names(colors) <- c("community","home", NA,  "other",NA,"social",NA,"work")
q <- ggplot(df.plot2, aes(x=tau, y=beta)) + 
  theme_bw(base_size = base.font.size) +
  geom_area(aes(fill=label)) +
  labs(x = expression(paste(tau, " (days)")),
       y = expression(paste(p, "(", tau, ")  (transmissions per day)")),
       title = NULL,  
       subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, 0.5), expand = F) +
  scale_fill_manual(values = colors[c("home", "work",  "social","community","other")],
                    labels = c(
                      bquote(paste('Dwelling')),
                      bquote(paste('Work')),
                      bquote(paste('Cospace-time')),
                      bquote(paste('Community')),
                      bquote(paste('Other')))) +
  scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks)

q
plot.name <- paste0("supple_figure_ptau.pdf")
ggsave(plot.name, q, height = plot.height, width = 11.75)


###export ptau for contact location of dwelling for calculations of secondary cases
write.table(df.p.tau.loc.c1,"ptau_by_time.csv",row.names = FALSE,col.names = TRUE,sep = ",")
###export contributions of transmission under five locations
frac<-as.data.frame(matrix(data = NA,nrow = 5,ncol = 2))
colnames(frac)<-c("patt","contribution")
frac[,1]<-c("Dwelling", "Work",  "Cospace-time","Community","Other")
frac[,2]<-c(frac_c1,frac_c2,frac_c3,frac_c4,frac_c5)
write.table(frac,"frac_of_R_across_all_locs.csv",row.names = F,col.names = T,sep = ",")

