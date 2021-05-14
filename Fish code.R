# ========================================================================================================== #

library(chaoUtility)
library(readr)
library(ape)
library(dplyr)
library(ggplot2)
library(Rcpp)
library(reshape2)
library(ggpubr)
library(gg.gap)


# ========================================================================================================== #

library(devtools)
install_github("AnneChao/iNEXT3D")  # Press 'Enter' to skip number selection
library(iNEXT3D)

source("figure-making.txt")

# ========================================================================================================== #
# Part 1 : Abundance-based yearly analysis for Figure 1, 2.

Abun <- read.csv("Fish Abundance data.csv", row.names = 1, header= TRUE)
tree <- read.tree("Fish PhyTree.txt")
traits <- read.csv("Fish Traits.csv", row.names = 1, header= TRUE)


Cmax <- apply(Abun, 2, function(x) iNEXT3D:::Chat.Ind(x, 2*sum(x))) %>% min %>% round(., 4)
Cmin <- apply(Abun, 2, function(x) iNEXT3D:::Chat.Ind(x, sum(x))) %>% min %>% round(., 4)


# ========================================================================================================== #
# Figure 1 - Taxonomic diversity

TD_est <- estimate3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', base = 'coverage', level = c(Cmin, Cmax), nboot = 0)
TD_obs <- Obs3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
TD_asy <- Asy3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0)


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qD, goalSC), 
                TD_obs %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed') , 
                TD_asy %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic') )
fig_1_or_3(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# Figure 1 - Phylogenetic diversity

PD_est <- estimate3D(data = Abun, diversity = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD', 
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, reftime = 1)
PD_obs <- Obs3D(data = Abun, diversity = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD',
                q = c(0, 1, 2), nboot = 0, reftime = 1) 
PD_asy <- Asy3D(data = Abun, diversity = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD',
                q = c(0, 1, 2), nboot = 0, reftime = 1)


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, goalSC), 
                PD_obs %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic') )
fig_1_or_3(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# Figure 1 - Functional diversity

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[,i] <- factor(traits[,i], levels = unique(traits[,i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Abun, diversity = 'FD', distM = distM, q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, datatype = 'abundance')
FD_obs <- Obs3D(data = Abun, diversity = 'FD', distM = distM, q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
FD_asy <- Asy3D(data = Abun, diversity = 'FD', distM = distM, datatype = 'abundance', q = c(0, 1, 2), nboot = 0)


out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD = qAUC, goalSC), 
                FD_obs %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed') , 
                FD_asy %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic'))
fig_1_or_3(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# Figure 2

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 0)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 1)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 2)


# ========================================================================================================== #
# Part 2 : Incidence-based three years analysis for Figure 3, 4.

Inci_freq <- read.csv("Fish Incidence frequency data.csv", row.names = 1, header= TRUE)
Inci_raw <- read.csv("Fish Incidence raw data.csv", row.names = 1, header= TRUE)
nT <- as.matrix(Inci_freq)[1,]
tree <- read.tree("Fish PhyTree.txt")
traits <- read.csv("Fish Traits.csv", row.names = 1, header= TRUE)


Cmax <- apply(Inci_freq, 2, function(x) iNEXT3D:::Chat.Sam(x, 2*x[1])) %>% min %>% round(., 4)
Cmin <- apply(Inci_freq, 2, function(x) iNEXT3D:::Chat.Sam(x, x[1])) %>% min %>% round(., 4)


# ========================================================================================================== #
# Figure 3 - Taxonomic diversity

TD_est <- estimate3D(data = Inci_freq, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', base = 'coverage', 
                     level = c(Cmin, Cmax), nboot = 0)
TD_obs <- Obs3D(data = Inci_freq, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', nboot = 0)
TD_asy <- Asy3D(data = Inci_freq, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', nboot = 0)


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qD, goalSC), 
                TD_obs %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed'), 
                TD_asy %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic'))
fig_1_or_3(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# Figure 3 - Phylogenetic diversity

PD_est <- estimate3D(data = Inci_raw, diversity = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, reftime = 1)
PD_obs <- Obs3D(data = Inci_raw, diversity = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                q = c(0, 1, 2), nboot = 0, reftime = 1)
PD_asy <- Asy3D(data = Inci_raw, diversity = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                q = c(0, 1, 2), nboot = 0, reftime = 1)


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, goalSC), 
                PD_obs %>% select(Assemblage, Order.q, qPD ) %>% mutate(goalSC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic'))
fig_1_or_3(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# Figure 3 - Functional diversity

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[, i] <- factor(traits[,i], levels = unique(traits[, i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Inci_freq, diversity = 'FD', distM = distM, datatype = 'incidence_freq',
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0)
FD_obs <- Obs3D(data = Inci_freq, diversity = 'FD', distM = distM, datatype = 'incidence_freq', 
                q = c(0, 1, 2), nboot = 0)
FD_asy <- Asy3D(data = Inci_freq, diversity = 'FD', distM = distM, datatype = 'incidence_freq',
                q = c(0, 1, 2), nboot = 0)

out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD = qAUC, goalSC), 
                FD_obs %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed'), 
                FD_asy %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic'))
fig_1_or_3(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# Figure 4

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 0)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 1)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 2)


