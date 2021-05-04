# ========================================================================================================== #
# 
# R code for computing three diversity of data fish in the paper:
# Chao, A., Henderson, P. A., (2021). 
# Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT-3D standardization.
# 
# This code includes two parts:
# (1) Abundance-based yearly analysis.
# (2) Incidence-based three-year analysis.
# 
# Please user install and library following relative packages first.
# 
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
# 
# Before plotting diversity, we need to introduce the most important package : 'iNEXT3D'.
# It is a package which is used to calculate Taxonomic Diversity, Phylogenetic Diversity, and Functional Diversity. 
# The main function in this package are 'iNEXT3D', 'estimate3D', 'Asy3D', and 'Obs3D'.
# 'iNEXT3D' is used to calculate interpolation and extrapolation diversity under specified sample size.
# 'estimate3D' is used to calculate interpolation and extrapolation diversity under specified sample coverage.
# 'Asy3D' is used to calculate asymptotic diversity. 'Obs3D' is used to calculate empirical diversity.
# In this paper, we will express four standardized diversity: 
# 'Asymptotic diversity', 'empirical diversity', and 'estimate3D' which can be calculated diversity under Cmin and Cmax. 
# If you want to get more details, please key in '?iNEXTeD', '?estimate3D', '?Asy3D', '?Obs3D' to get more information.
# Finally, we can compare different time series pattern under any levels.
# 
# ========================================================================================================== #

library(devtools)
install_github("AnneChao/iNEXT3D")  # Press 'Enter' to skip number selection
library(iNEXT3D)


# ========================================================================================================== #
# 
# Except packages, we also need to run following two plotting functions :
# (1) Compare diversity versus time series for all order q under a specified diversity class. (For figure 1, 3)
# (2) Compare diversity versus time series for three diversity class under a specified order q. (For figure 2, 4)
# 
# ========================================================================================================== #

#' @param fig_13 is used to plot figure 1, figure 3 in the paper
#' @param output is the result come from iNEXT3D.
#' @param y_label is the diversity type.
#' 
fig_13 <- function(output, y_label) {
  goalSC = unique(output$goalSC)[1:2] %>% as.numeric
  output <- as_tibble(output)
  output$goalSC <- factor(output$goalSC, levels = c(goalSC[1], 'Observed', goalSC[2], 'Asymptotic'))
  poly_ord <- c(1, 4)
  q = unique(output$Order.q)
  fit_info <- matrix(NA, nrow = length(poly_ord)*length(q), ncol = 2*nlevels(output$goalSC))
  names(output)[3] <- 'qD'
  for(i in 1:length(q)){
    for(j in 1:nlevels(output$goalSC)){
      myout_ <- output %>% filter(goalSC == levels(output$goalSC)[j], Order.q == q[i]) %>% 
        mutate(year = as.numeric(substr(Assemblage, 3, 8))-1980)
      for(k in 1:length(poly_ord)){
        tmpp <- lm(formula = qD ~ poly(year,poly_ord[k]), data = myout_) %>% summary
        anov <- tmpp$coefficients[nrow(tmpp$coefficients), ncol(tmpp$coefficients)] %>%
          round(.,4)
        rsq_adj <- tmpp$adj.r.squared %>% round(.,4)
        fit_info[length(poly_ord)*i-(k-1), c(2*j-1,2*j)] <- c(anov, rsq_adj)
      }
    }
  }
  colnames(fit_info) <- rep(levels(output$goalSC), each = 2)
  rownames(fit_info) <- rep(paste0('q = ', q), each = length(poly_ord))
  fit_info2 <- fit_info[, rev(1:ncol(fit_info))]
  anovs <- fit_info[, seq(1, ncol(fit_info), 2)]
  anovs <- cbind(Order.q = rep(q,each = length(poly_ord)), poly_ord = rep(rev(poly_ord), length(q)), anovs) %>% 
    as_tibble(.) %>% 
    melt(.,id.vars = c('Order.q', 'poly_ord'), variable.name = 'goalSC',value.name = 'pvalue') %>% as_tibble() %>% 
    mutate(sig = as.numeric(pvalue < 0.05)) %>% select(-pvalue)
  
  output <- output %>% mutate(year = as.numeric(substr(Assemblage, 3, nchar(Assemblage))) - 1980)
  output <- output %>% arrange(goalSC, Order.q)
  output <- output %>% group_by(goalSC, Order.q) %>% 
    do(lm(formula = qD ~ poly(year,4), data = . ) %>% 
         predict %>% tibble(fitD4 = .)) %>% 
    bind_cols(output) %>% ungroup %>% select(year, Order.q, qD, goalSC, fitD4)
  output <- output %>% group_by(goalSC, Order.q) %>% 
    do(lm(formula = qD ~ poly(year, 1), data = . ) %>% 
         predict %>% tibble(fitD1 = .)) %>% 
    bind_cols(output) %>% ungroup %>% select(year, Order.q, qD, goalSC, fitD4, fitD1)
  
  output <- melt(output, id.vars = c('year','Order.q','goalSC'), variable.name = 'type', value.name = 'qD') %>% 
    as_tibble() 
  output <- output %>% mutate(poly_ord = ifelse(type == 'qD', 0, ifelse(type == 'fitD4', 4, 1)))
  output <- left_join(x = output, y = anovs, by = c('Order.q', 'goalSC', 'poly_ord')) %>%
    select(-poly_ord)
  output$sig[is.na(output$sig)] <- 0
  output$goalSC <- as.character(output$goalSC)
  
  output$goalSC[output$goalSC != 'Asymptotic' & output$goalSC != 'Observed'] <- 
    paste0('Coverage = ', output$goalSC[output$goalSC != 'Asymptotic' & output$goalSC != 'Observed'] )
  output$type <- as.character(output$type)
  
  ##
  pics <- list()
  goalSC_title <- rev(unique(output$goalSC))
  goalSC_title_tmp <- goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"]
  goalSC_title_tmp <- goalSC_title_tmp %>% gsub(pattern = 'Coverage = ', replacement = '', x = .) %>% 
    as.numeric(.)*100 
  goalSC_title_tmp <- 
    sapply(goalSC_title_tmp, function(x) {
      ifelse(round(x)-x == 0,substr(as.character(x), 1, 2), substr(as.character(x), 1, nchar(x)))
    }) %>% paste0('Coverage = ',.,'%')
  goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"] <- goalSC_title_tmp
  output <- output %>% filter(goalSC %in% rev(unique(output$goalSC)))
  output$goalSC <- factor(output$goalSC, levels = rev(unique(output$goalSC)))
  output$type <- factor(output$type, levels = unique(output$type))
  maxy <- max(output$qD)
  miny <- min(output$qD)
  
  if (y_label == 'Taxonomic diversity') {
    n_break = 5
  } else if (y_label == 'Phylogenetic diversity') {
    n_break = 6
  } else {n_break = 7}
  
  for (i in 1:nlevels(output$goalSC)){
    tmp <- output %>% filter(goalSC == levels(output$goalSC)[i])
    if(length(unique(tmp$sig)) == 1){
      tmp <- rbind(tmp, tmp[1,])
      tmp[nrow(tmp), "qD"] <- NA
      tmp[nrow(tmp), "sig"] <- 1-unlist(tmp[nrow(tmp), "sig"])
    }
    tmp$goalSC <- as.character(tmp$goalSC)
    tmp$sig <- factor(tmp$sig, levels = c(0, 1))
    tmp$Order.q <- factor(tmp$Order.q, levels = q)
    
    #tmp <- tmp %>% group_by(type)
    pp <- ggplot(data = tmp) + theme_bw() +
      geom_line(aes(x = year + 1980, y = qD, size = sig, colour = Order.q, alpha = type, linetype = type)) +
      coord_cartesian(ylim = c(miny, maxy)) +
      scale_linetype_manual(values = c(1, 1, 2), guide = FALSE) +
      scale_color_manual(values = c('#1F78B4', '#E7298A', '#1B9E77'), name = "Order.q") +
      scale_alpha_manual(values = c(0.4, 1, 1), guide = FALSE) +
      scale_size_manual(values = c(1, 1.6), guide = FALSE) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = n_break)) +
      theme(legend.position = 'bottom', legend.text = element_text(size = 18),
            legend.title = element_text(size = 18),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 9),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 17)) +
      ggtitle(goalSC_title[i]) + 
      guides(color = guide_legend(override.aes = list(size = 2)))
    pics[[i]] <- pp
  }
  ans <- ggarrange(plotlist = pics, ncol = 4, nrow = 1, common.legend = TRUE, legend = 'bottom')
  ans <- annotate_figure(ans, left = text_grob(y_label, rot = 90, size = 18, hjust = 0.5))
  ans
}

#' @param fig_24 is used to plot figure 2, figure 3 in the paper
#' @param TD.output output from taxonomic class.
#' @param PD.output output from phylogenetic class.
#' @param FD.output output from functional class.
#' @param q is the specified order q.
#' 
fig_24 <- function(TD.output, PD.output, FD.output, q) {
  find_ord <- function(myout) {
    goalSC = unique(myout$goalSC)[1:2] %>% as.numeric
    q = unique(myout$Order.q)
    myout <- as_tibble(myout)
    myout$goalSC <- factor(myout$goalSC, levels = c(goalSC[1], 'Observed', goalSC[2], 'Asymptotic'))
    poly_ord <- c(1, 4)
    fit_info <- matrix(NA, nrow = length(poly_ord)*length(q), ncol = 2*nlevels(myout$goalSC))
    names(myout)[3] <- 'qD'
    for(i in 1:length(q)){
      for(j in 1:nlevels(myout$goalSC)){
        myout_ <- myout %>% filter(goalSC == levels(myout$goalSC)[j], Order.q == q[i]) %>% 
          mutate(year = as.numeric(substr(Assemblage, 3, 8))-1980)
        for(k in 1:length(poly_ord)){
          tmpp <- lm(formula = qD ~ poly(year,poly_ord[k]), data = myout_) %>% summary
          anov <- tmpp$coefficients[nrow(tmpp$coefficients), ncol(tmpp$coefficients)] %>%
            round(.,4)
          rsq_adj <- tmpp$adj.r.squared %>% round(.,4)
          fit_info[length(poly_ord)*i-(k-1), c(2*j-1,2*j)] <- c(anov, rsq_adj)
        }
      }
    }
    colnames(fit_info) <- rep(levels(myout$goalSC), each = 2)
    rownames(fit_info) <- rep(paste0('q = ', q), each = length(poly_ord))
    fit_info2 <- fit_info[, rev(1:ncol(fit_info))]
    anovs <- fit_info[, seq(1, ncol(fit_info), 2)]
    anovs <- cbind(Order.q = rep(q,each = length(poly_ord)), poly_ord = rep(rev(poly_ord), length(q)), anovs) %>% 
      as_tibble(.) %>% 
      melt(.,id.vars = c('Order.q', 'poly_ord'), variable.name = 'goalSC',value.name = 'pvalue') %>% as_tibble() %>% 
      mutate(sig = as.numeric(pvalue < 0.05)) %>% select(-pvalue)
    
    myout <- myout %>% mutate(year = as.numeric(substr(Assemblage, 3, nchar(Assemblage))) - 1980)
    myout <- myout %>% arrange(goalSC, Order.q)
    myout <- myout %>% group_by(goalSC, Order.q) %>% 
      do(lm(formula = qD ~ poly(year,4), data = . ) %>% 
           predict %>% tibble(fitD4 = .)) %>% 
      bind_cols(myout) %>% ungroup %>% select(year, Order.q, qD, goalSC, fitD4)
    myout <- myout %>% group_by(goalSC, Order.q) %>% 
      do(lm(formula = qD ~ poly(year, 1), data = . ) %>% 
           predict %>% tibble(fitD1 = .)) %>% 
      bind_cols(myout) %>% ungroup %>% select(year, Order.q, qD, goalSC, fitD4, fitD1)
    
    myout <- melt(myout, id.vars = c('year','Order.q','goalSC'), variable.name = 'type', value.name = 'qD') %>% 
      as_tibble() 
    myout <- myout %>% mutate(poly_ord = ifelse(type == 'qD', 0, ifelse(type == 'fitD4', 4, 1)))
    myout <- left_join(x = myout, y = anovs, by = c('Order.q', 'goalSC', 'poly_ord')) %>%
      select(-poly_ord)
    myout$sig[is.na(myout$sig)] <- 0
    myout$goalSC <- as.character(myout$goalSC)
    
    myout$goalSC[myout$goalSC != 'Asymptotic' & myout$goalSC != 'Observed'] <- 
      paste0('Coverage = ', myout$goalSC[myout$goalSC != 'Asymptotic' & myout$goalSC != 'Observed'] )
    myout$type <- as.character(myout$type)
    myout
  }
  
  out_td <- find_ord(myout = TD.output) %>% mutate(type2 = 'Taxonomic')
  out_pd <- find_ord(myout = PD.output) %>% mutate(type2 = 'Phylogenetic')
  out_fd <- find_ord(myout = FD.output) %>% mutate(type2 = 'Functional')
  output <- rbind(out_td, out_pd, out_fd)
  
  pics <- list()
  goalSC_title <- rev(unique(output$goalSC))
  goalSC_title_tmp <- goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"]
  goalSC_title_tmp <- goalSC_title_tmp %>% gsub(pattern = 'Coverage = ', replacement = '', x = .) %>% 
    as.numeric(.)*100 
  goalSC_title_tmp <- 
    sapply(goalSC_title_tmp, function(x) {
      ifelse(round(x)-x == 0, substr(as.character(x), 1, 2), substr(as.character(x), 1, nchar(x)))
    }) %>% paste0('Coverage = ', ., '%')
  goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"] <- goalSC_title_tmp
  output <- output %>% filter(goalSC %in% rev(unique(output$goalSC)), Order.q == q)
  output$goalSC <- factor(output$goalSC, levels = rev(unique(output$goalSC)))
  output$type <- factor(output$type, levels = unique(output$type))
  
  floor_dec <- function(x, dig = 1) round(x - 5*10^(- dig - 1), dig)
  ceiling_dec <- function(x, dig = 1) round(x + 5*10^(- dig - 1), dig)
  
  maxy <- ceiling(max(output$qD))
  
  if (q == 0) {
    miny <- floor(min(output$qD))
  } else {
    miny <- floor_dec(min(output$qD), 1)
  }
  
  if (q == 0) {
    y_max_low <- ceiling(max(output$qD[output$type2 != 'Taxonomic'], na.rm = T))
    y_min_high <- floor(min(output$qD[output$type2 == 'Taxonomic'], na.rm = T))
  }
  
  for (i in 1:nlevels(output$goalSC)) {
    tmp <- output %>% filter(goalSC == levels(output$goalSC)[i])
    if (length(unique(tmp$sig)) == 1){
      tmp <- rbind(tmp, tmp[1,])
      tmp[nrow(tmp), "qD"] <- NA
      tmp[nrow(tmp), "sig"] <- 1 - unlist(tmp[nrow(tmp), "sig"])
    }
    
    tmp$goalSC <- as.character(tmp$goalSC)
    tmp$sig <- factor(tmp$sig, levels = c(0, 1))
    tmp$type2 <- factor(tmp$type2, levels = unique(tmp$type2))
    
    pp <- ggplot(data = tmp) + theme_bw() +
      geom_line(aes(x = year + 1980, y = qD, size = sig, colour = type2, alpha = type, linetype = type))+
      coord_cartesian(ylim = c(miny, maxy))+
      scale_linetype_manual(values = c(1, 1, 2), guide = FALSE)+
      scale_color_manual(values = c('#A700D5', 'black', '#D55E00'))+
      scale_alpha_manual(values = c(0.4, 1, 1), guide = FALSE)+
      scale_size_manual(values = c(1, 1.8), guide = FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
      theme(legend.position = 'bottom', legend.text = element_text(size = 18),
            legend.title = element_blank(),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, size = 17)) +
      ggtitle(goalSC_title[i]) + ylab(NULL) + xlab(NULL) + 
      guides(color = guide_legend(override.aes = list(size = 2)))
    
    if(q == 0){
      tick_width <- c(round((y_max_low - miny)/2.5), round((maxy - y_min_high)/2.5))
      pp <- gg.gap(plot = pp, ylim = c(miny, maxy),
                   tick_width = tick_width,
                   segments = c(y_max_low, y_min_high))+
        theme(plot.margin = unit(c(1, 0, 1, 0), "lines"), legend.position = 'bottom')
    }
    pics[[i]] <- pp
  }
  ans <- ggarrange(plotlist = pics, ncol = 4, nrow = 1,
                   common.legend = TRUE, legend = 'bottom')
  ans <- annotate_figure(ans, left = text_grob('Diversity', rot = 90, size = 17, hjust = 0.5))
  ans
}


# ========================================================================================================== #
# 
# Part 1 : Abundance-based yearly analysis for Figure 1, 2.
# 
# For abundance yearly analysis, we need to load abundance data, phylogenetic tree, and traits.
# Abundance data is composed of species ~ goup.
# Each row means a species and each column means a group.
# Each numeric means abundances frequency of each species in a assemblage. 
# 
# Phylogenetic tree is a 'phylo' class data which we can know the evolution time of each species.
# Traits is a data which record the average traits of each species.
# 
# And then we need to compute Cmin and Cmax.
# Cmin : the smallest observed coverage among all assemblages.
# Cmax : the smallest coverage after extrapolating double sample size among all assemblages.
# 
# ========================================================================================================== #

Abun <- read.csv("Fish Abundance data.csv", row.names = 1, header= TRUE)
tree <- read.tree("Fish PhyTree.txt")
traits <- read.csv("Fish Traits.csv", row.names = 1, header= TRUE)


Cmax <- apply(Abun, 2, function(x) iNEXT3D:::Chat.Ind(x, 2*sum(x))) %>% min %>% round(., 4)
Cmin <- apply(Abun, 2, function(x) iNEXT3D:::Chat.Ind(x, sum(x))) %>% min %>% round(., 4)


# ========================================================================================================== #
# 
# Figure 1 - Taxonomic diversity
# Here we only use abundance data. 
# First, calculate taxonomic diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

TD_est <- estimate3D(data = Abun, class = 'TD', q = c(0, 1, 2), datatype = 'abundance', base = 'coverage', level = c(Cmin, Cmax), nboot = 0)
TD_obs <- Obs3D(data = Abun, class = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
TD_asy <- Asy3D(data = Abun, class = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0)


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qD, goalSC), 
                TD_obs %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed') , 
                TD_asy %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic') )
fig_13(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# 
# Figure 1 - Phylogenetic diversity
# For phylogenetic diversity, we will use abundance data and phylogenetic tree.
# First, calculate phylogenetic diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

PD_est <- estimate3D(data = Abun, class = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD', 
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, reftime = 1)
PD_obs <- Obs3D(data = Abun, class = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD',
                q = c(0, 1, 2), nboot = 0, reftime = 1) 
PD_asy <- Asy3D(data = Abun, class = 'PD', tree = tree, datatype = 'abundance', PDtype = 'meanPD',
                q = c(0, 1, 2), nboot = 0, reftime = 1)


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, goalSC), 
                PD_obs %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic') )
fig_13(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# 
# Figure 1 - Functional diversity
# For functional diversity, we will use abundance data and transform traits into distance matrix by metric 'gower' distance.
# Before transforming, we should transform character variable into factor variable. 
# The value of distance matrix should between 'zero' and 'one'.
# First, calculate functional diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[,i] <- factor(traits[,i], levels = unique(traits[,i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Abun, class = 'AUC', distM = distM, q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, datatype = 'abundance')
FD_obs <- Obs3D(data = Abun, class = 'AUC', distM = distM, q = c(0, 1, 2), datatype = 'abundance', nboot = 0)
FD_asy <- Asy3D(data = Abun, class = 'AUC', distM = distM, datatype = 'abundance', q = c(0, 1, 2), nboot = 0)


out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD = qAUC, goalSC), 
                FD_obs %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed') , 
                FD_asy %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic'))
fig_13(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# 
# Figure 2
# Here we also use the previous result to plot diversity with three class diversity groups under a specified order q. 
# 
# ========================================================================================================== #

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 0)

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 1)

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 2)


# ========================================================================================================== #
# 
# Part 2 : Incidence-based three years analysis for Figure 3, 4.
# 
# For incidence three years analysis, we need to load incidence data, phylogenetic tree, and traits.
# Incidence data can be classified into 'incidence-raw data' and 'incidence-frequency data'.
# 'Incidence-raw data' can be got by 'incidence-frequency data' with specified sampling units.
# 
# In 'incidence-frequency data', it is composed of species ~ group.
# Each row means a species besides first row which means sampling units and each column means a group.
# Except for first row, each numeric means how many sampling units a species is found in a assemblage.
# Thus under a certain assemblage, the incidences should not be larger than total sampling units.
# 
# In 'incidence-raw data', it is composed of species ~ sampling unit.
# Each numeric means whether the species is found in this sampling unit.
# Thus, the value in 'incidence-raw data' is only '0' (not found) or '1' (found).
# User should combine each group by column to get a big matrix. And key in how many sampling units in each group ('nT').
# 
# Here we will use 'incidence-frequency data' to compute 'Taxonomic diversity' and 'Functional diversity',
# and also use 'incidence-raw data' to compute 'Phylogenetic diversity'.
# 
# Phylogenetic tree is a 'phylo' class data which we can know the evolution time of each species.
# Traits is a data which record the average traits of each species.
# 
# And then we need to compute Cmin and Cmax.
# Cmin : the smallest observed coverage among all assemblages.
# Cmax : the smallest coverage after extrapolating double sample size among all assemblages.
# 
# ========================================================================================================== #

Inci_freq <- read.csv("Fish Incidence frequency data.csv", row.names = 1, header= TRUE)
Inci_raw <- read.csv("Fish Incidence raw data.csv", row.names = 1, header= TRUE)
nT <- as.matrix(Inci_freq)[1,]
tree <- read.tree("Fish PhyTree.txt")
traits <- read.csv("Fish Traits.csv", row.names = 1, header= TRUE)


Cmax <- apply(Inci_freq, 2, function(x) iNEXT3D:::Chat.Sam(x, 2*x[1])) %>% min %>% round(., 4)
Cmin <- apply(Inci_freq, 2, function(x) iNEXT3D:::Chat.Sam(x, x[1])) %>% min %>% round(., 4)


# ========================================================================================================== #
# 
# Figure 3 - Taxonomic diversity
# Here we only use incidence-frequency data. 
# First, calculate taxonomic diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

TD_est <- estimate3D(data = Inci_freq, class = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', base = 'coverage', 
                     level = c(Cmin, Cmax), nboot = 0)
TD_obs <- Obs3D(data = Inci_freq, class = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', nboot = 0)
TD_asy <- Asy3D(data = Inci_freq, class = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', nboot = 0)


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qD, goalSC), 
                TD_obs %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed'), 
                TD_asy %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic'))
fig_13(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# 
# Figure 3 - Phylogenetic diversity
# For phylogenetic diversity, we will use incidence-raw data and phylogenetic tree.
# First, calculate phylogenetic diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

PD_est <- estimate3D(data = Inci_raw, class = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0, reftime = 1)
PD_obs <- Obs3D(data = Inci_raw, class = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                q = c(0, 1, 2), nboot = 0, reftime = 1)
PD_asy <- Asy3D(data = Inci_raw, class = 'PD', tree = tree, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                q = c(0, 1, 2), nboot = 0, reftime = 1)


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, goalSC), 
                PD_obs %>% select(Assemblage, Order.q, qPD ) %>% mutate(goalSC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic'))
fig_13(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# 
# Figure 3 - Functional diversity
# For functional diversity, we will use incidence-frequency data and transform traits into distance matrix by metric 'gower' distance.
# Before transforming, we should transform character variable into factor variable. 
# The value of distance matrix should between 'zero' and 'one'.
# First, calculate functional diversity under Cmin and Cmax, asymptotic diversity and empirical diversity by 'iNEXT3D'.
# Then, we can plot it with group 'q = 0', 'q = 1', 'q = 2' under each standardized coverage.
# 
# ========================================================================================================== #

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[, i] <- factor(traits[,i], levels = unique(traits[, i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Inci_freq, class = 'AUC', distM = distM, datatype = 'incidence_freq',
                     q = c(0, 1, 2), level = c(Cmin, Cmax), nboot = 0)
FD_obs <- Obs3D(data = Inci_freq, class = 'AUC', distM = distM, datatype = 'incidence_freq', 
                q = c(0, 1, 2), nboot = 0)
FD_asy <- Asy3D(data = Inci_freq, class = 'AUC', distM = distM, datatype = 'incidence_freq',
                q = c(0, 1, 2), nboot = 0)

out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD = qAUC, goalSC), 
                FD_obs %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed'), 
                FD_asy %>% select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic'))
fig_13(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# 
# Figure 4
# Here we also use the previous result to plot diversity with three class diversity groups under a specified order q. 
# 
# ========================================================================================================== #

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 0)

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 1)

fig_24(TD = out_TD, PD = out_PD, FD = out_FD, q = 2)


