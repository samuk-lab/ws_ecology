######################################################################
# Analysis of isotopic abundance data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("ggplot2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("SIBER")

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

setwd("analysis_morphology")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
iso_df <- read.csv("data/isotope_data.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

iso_df  <- left_join(iso_df, meta_df, by = "id")

iso_df <- iso_df %>%
  mutate(group = ifelse(cluster == "wht", "white", "common")) %>%
  mutate(pop = gsub("[^A-Z]*", "", as.character(id))) %>%
  mutate(cn.ratio = C.amount / N.amount)

# reformat for SIBER

#LN and AL are monospecies sites, so removing

iso_df_siber <- iso_df %>%
  filter(!(pop %in% c("LN", "AL"))) %>%
  filter(!is.na(group)) %>%
  select(d13C, d15N, group, pop) 

names(iso_df_siber) <- c("iso1", "iso2", "group", "community")

data("demo.siber.data")

siber.example <- createSiberObject(iso_df_siber)

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# models

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(siber.example, parms, priors)


