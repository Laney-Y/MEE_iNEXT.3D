library(chaoUtility)
library(readr)
library(ape)
library(dplyr)
library(ggplot2)
library(Rcpp)
library(reshape2)
library(ggpubr)
library(gg.gap)

library(devtools)
# install_github("AnneChao/iNEXT3D")
library(iNEXT3D)

######## functions -------------------------------------------------------------------
## Data cleaning : following functions 'abun', 'inci_freq', and 'inci_raw', transform monthly abundance data into abundance-based, incidence-frequency-based, and incidence-raw-based data.
# Arguments : 
# 'dataset' is monthly data profile.
# 'frwidth' is the width range of year in each group.
abun <- function(dataset, frwidth = 1) {
  dataset[is.na(dataset)] <- 0
  dataset <- dataset[which(rowSums(dataset[,-1]) > 0),]
  dataset <- dataset %>% mutate(year = as.numeric(substr(dataset$`Date of Sample`, start = 7, stop = 10)))
  year <- dataset$year
  tmp <- dataset[, -c(1, ncol(dataset))]
  dataset <- cbind(year = dataset$year, tmp)
  colnames(dataset)[-c(1)] <- gsub(x = colnames(dataset)[-c(1)], pattern = " ", replacement = "_")
  
  dataset <- dataset %>% group_by(year) %>% summarise_at(colnames(dataset)[-1], sum) %>% t()
  colnames(dataset) <- paste0('yr', dataset[1,]); dataset <- dataset[-1,]
  
  dataset <- dataset[, -which(colnames(dataset) == 'yr1986')]
  
  old_sp <- read_csv("Fish Abundance data.csv") %>% .[,1] %>% unlist()
  rownames(dataset)[which(rownames(dataset) == "Gaidropsaurus_vulgaris")] <- 'Gaidropsarus_vulgaris'
  rownames(dataset)[which(rownames(dataset) == "Lampetra_fluviatalis")] <- 'Lampetra_fluviatilis'
  rownames(dataset)[which(rownames(dataset) == "Chelon_auritus")] <- 'Chelon_auratus'
  rownames(dataset)[which(rownames(dataset) == "Scophalmus_maximus")] <- 'Scophthalmus_maximus'
  dataset <- dataset[-which(rownames(dataset)=='Perca_fluviatilis'), ]
  print(sum(rownames(dataset) %in% old_sp) == length(old_sp))
  
  names_drop <- intersect(rowSums(dataset)[rowSums(dataset) > 1000] %>% names(),
                          rowSums(dataset)[rowSums(dataset) > 1000] %>% names())
  
  if(frwidth > 1){
    yr_mean = sapply(1:(ncol(dataset)-frwidth+1), function(i){
      colnames(dataset[, seq(i, i+frwidth-1, 1)]) %>% 
        substr(., start = 3, stop = nchar(.)) %>% as.numeric() %>% mean %>% round(., 1)
    }) %>% paste0('yr',.)
    
    dataset <- sapply(1:(ncol(dataset)-frwidth+1), function(i){
      rowSums(dataset[, seq(i, i+frwidth-1, 1)])
    })
    
    colnames(dataset) <- yr_mean
  }
  dataset_abun <- dataset[match(names_drop, rownames(dataset)),]
  dataset_rare <- dataset[-match(names_drop, rownames(dataset)),]
  return(list(abun = dataset_abun, rare = dataset_rare))
}
inci_freq <- function(dataset, frwidth = 1) {
  dataset[is.na(dataset)] <- 0
  dataset <- dataset[which(rowSums(dataset[,-1]) > 0),]
  dataset <- dataset %>% mutate(year = as.numeric(substr(dataset$`Date of Sample`, start = 7, stop = 10)))
  year <- dataset$year
  tmp <- dataset[, -c(1,ncol(dataset))];tmp[tmp>1] <- 1
  dataset <- cbind(year = dataset$year, nT = 1, tmp)
  colnames(dataset)[-c(1,2)] <- gsub(x = colnames(dataset)[-c(1,2)], pattern = " ", replacement = "_")
  
  dataset <- dataset %>% group_by(year) %>% summarise_at(colnames(dataset)[-1], sum) %>% t()
  colnames(dataset) <- paste0('yr', dataset[1,]); dataset <- dataset[-1,]
  
  dataset <- dataset[, -which(colnames(dataset) == 'yr1986')]
  
  old_sp <- read_csv("Fish Abundance data.csv") %>% .[,1] %>% unlist()
  rownames(dataset)[which(rownames(dataset) == "Gaidropsaurus_vulgaris")] <- 'Gaidropsarus_vulgaris'
  rownames(dataset)[which(rownames(dataset) == "Lampetra_fluviatalis")] <- 'Lampetra_fluviatilis'
  rownames(dataset)[which(rownames(dataset) == "Chelon_auritus")] <- 'Chelon_auratus'
  rownames(dataset)[which(rownames(dataset) == "Scophalmus_maximus")] <- 'Scophthalmus_maximus'
  dataset <- dataset[-which(rownames(dataset)=='Perca_fluviatilis'),]
  print(sum(rownames(dataset) %in% old_sp) == length(old_sp))
  
  if(frwidth > 1){
    yr_mean = sapply(1:(ncol(dataset)-frwidth+1), function(i) {
      colnames(dataset[, seq(i,i+frwidth-1,1)]) %>%
        substr(., start = 3, stop = nchar(.)) %>% as.numeric() %>% mean %>% round(., 1)
    }) %>% paste0('yr',.)
    
    dataset <- sapply(1:(ncol(dataset)-frwidth+1), function(i){
      rowSums(dataset[, seq(i, i + frwidth - 1,1)])
    })
    
    colnames(dataset) <- yr_mean
  }
  dataset
}
inci_raw <- function(dataset, frwidth = 1) {
  dataset[is.na(dataset)] <- 0
  dataset <- dataset[which(rowSums(dataset[,-1]) > 0),]
  dataset <- dataset %>% mutate(year = as.numeric(substr(dataset$`Date of Sample`, start = 7, stop = 10)))
  year <- dataset$year
  tmp <- dataset[, -c(1,ncol(dataset))]; tmp[tmp>1] <- 1
  dataset <- cbind(year = dataset$year, nT = 1, tmp)
  colnames(dataset)[-c(1,2)] <- gsub(x = colnames(dataset)[-c(1,2)], pattern = " ",replacement = "_")
  
  dataset <- dataset[!(dataset$year == 1986),] %>% t()
  year = year[!(year == 1986)]
  
  old_sp <- read_csv("Fish Abundance data.csv") %>% .[,1] %>% unlist()
  rownames(dataset)[which(rownames(dataset) == "Gaidropsaurus_vulgaris")] <- 'Gaidropsarus_vulgaris'
  rownames(dataset)[which(rownames(dataset) == "Lampetra_fluviatalis")] <- 'Lampetra_fluviatilis'
  rownames(dataset)[which(rownames(dataset) == "Chelon_auritus")] <- 'Chelon_auratus'
  rownames(dataset)[which(rownames(dataset) == "Scophalmus_maximus")] <- 'Scophthalmus_maximus'
  dataset <- dataset[-which(rownames(dataset) == 'Perca_fluviatilis'), ]
  print(sum(rownames(dataset) %in% old_sp) == length(old_sp))
  
  yr_mean = sapply(1:(length(unique(year))-frwidth+1), function(i){
    unique(year)[seq(i, i+frwidth-1, 1)] %>% substr(., start = 1, stop = nchar(.)) %>% as.numeric() %>% mean %>% round(.,1)
  }) %>% paste0('yr',.)
  
  dataset <- lapply(1:(length(unique(year))-frwidth+1), function(i){
    dataset[-(1:2), dataset[1,] %in% unique(year)[seq(i, i+frwidth-1, 1)]]
  }) %>% do.call(cbind,.)
  
  dataset
}

##  Following function is used to check significance
find_ord <- function(myout) {
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

##  plot each diversity type versus to order q under along with sample coverage.
# Arguments :
# 'n_cols' is the number of figures,
# 'miny' and 'maxy' are the range of y scales.
# 'y_label' is the diversity type.
# 'nbreak' is the number of scales which will show on the figure.
# 'hjust_' is the horizontal adjustment of figures.
pic_qs <- function(myout, n_cols = 1:4, miny = NULL, maxy = NULL, y_label, nbreak = 5, hjust_ = 0.5) {
  pics <- list()
  goalSC_title <- rev(unique(myout$goalSC))[n_cols]
  goalSC_title_tmp <- goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"]
  goalSC_title_tmp <- goalSC_title_tmp %>% gsub(pattern = 'Coverage = ', replacement = '', x = .) %>% 
    as.numeric(.)*100 
  goalSC_title_tmp <- 
    sapply(goalSC_title_tmp, function(x) {
      ifelse(round(x)-x == 0,substr(as.character(x), 1, 2), substr(as.character(x), 1, nchar(x)))
    }) %>% paste0('Coverage = ',.,'%')
  goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"] <- goalSC_title_tmp
  myout <- myout %>% filter(goalSC %in% rev(unique(myout$goalSC))[n_cols])
  myout$goalSC <- factor(myout$goalSC, levels = rev(unique(myout$goalSC)))
  myout$type <- factor(myout$type, levels = unique(myout$type))
  if (is.null(miny)){
    maxy <- max(myout$qD)
    miny <- min(myout$qD)
  }
  for (i in 1:nlevels(myout$goalSC)){
    tmp <- myout %>% filter(goalSC == levels(myout$goalSC)[i])
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
      scale_y_continuous(breaks = scales::pretty_breaks(n = nbreak)) +
      theme(legend.position = 'bottom', legend.text = element_text(size = 18),
            legend.title = element_text(size = 18),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 11),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 17)) +
      ggtitle(goalSC_title[i])
    pics[[i]] <- pp
  }
  if (length(n_cols) <= 4){
    ans <- ggarrange(plotlist = pics, ncol = length(n_cols), nrow = 1,
                     common.legend = TRUE, legend = 'bottom')
  }else{
    ans <- ggarrange(plotlist = pics, ncol = 4, nrow = ceiling(length(n_cols)/4),
                     common.legend = TRUE, legend = 'bottom')
  }
  ans <- annotate_figure(ans, left = text_grob(y_label, rot = 90, size = 17, hjust = hjust_))
  ans
}

# plot three diversity type in the same plot under specified order q along with each sample coverage.
# Arguments :
# 'qfilter' is the specified order q.
# 'n_cols' is the number of figures,
# 'miny' and 'maxy' are the range of y scales.
# 'y_label' is the diversity type.
# 'nbreak' is the number of scales which will show on the figure.
# 'hjust_' is the horizontal adjustment of figures.
# 'seg' is the logical setting whether figure has segment or not.
# 'intg' is the logical setting whether y scales lower bound use minimum integer or not.
# 'max_2nd' is the logical setting whether y scales upper bound use second maximum integer or not.
pic_tpf <- function(myout, qfilter, n_cols = 1:4, y_label, nbreak = 5, hjust_ = 0.5, seg, intg, max_2nd = F) {
  pics <- list()
  goalSC_title <- rev(unique(myout$goalSC))[n_cols]
  goalSC_title_tmp <- goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"]
  goalSC_title_tmp <- goalSC_title_tmp %>% gsub(pattern = 'Coverage = ', replacement = '', x = .) %>% 
    as.numeric(.)*100 
  goalSC_title_tmp <- 
    sapply(goalSC_title_tmp, function(x) {
      ifelse(round(x)-x == 0, substr(as.character(x), 1, 2), substr(as.character(x), 1, nchar(x)))
    }) %>% paste0('Coverage = ', ., '%')
  goalSC_title[goalSC_title != "Observed" & goalSC_title != "Asymptotic"] <- goalSC_title_tmp
  myout <- myout %>% filter(goalSC %in% rev(unique(myout$goalSC))[n_cols], Order.q == qfilter)
  myout$goalSC <- factor(myout$goalSC, levels = rev(unique(myout$goalSC)))
  myout$type <- factor(myout$type, levels = unique(myout$type))
  
  floor_dec <- function(x, dig = 1) round(x - 5*10^(- dig - 1), dig)
  ceiling_dec <- function(x, dig = 1) round(x + 5*10^(- dig - 1), dig)
  
  if(max_2nd) {
    maxy <- ceiling(max(myout$qD[-which.max(myout$qD)]))
  } else {
    maxy <- ceiling(max(myout$qD))
  }
  
  if (intg) {
    miny <- floor(min(myout$qD))
  }else{
    miny <- floor_dec(min(myout$qD), 1)
  }
  
  if (seg) {
    if(intg){
      y_max_low <- ceiling(max(myout$qD[myout$type2 != 'Taxonomic'], na.rm = T))
    }else{
      y_max_low <- ceiling_dec(max(myout$qD[myout$type2 != 'Taxonomic'], na.rm = T), 1)
    }
    y_min_high <- floor(min(myout$qD[myout$type2 == 'Taxonomic'], na.rm = T))
  }
  for (i in 1:nlevels(myout$goalSC)) {
    tmp <- myout %>% filter(goalSC == levels(myout$goalSC)[i])
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
      scale_size_manual(values = c(1, 1.6), guide = FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = nbreak))+
      theme(legend.position = 'bottom', legend.text = element_text(size = 18),
            legend.title = element_blank(),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, size = 17))+
      ggtitle(goalSC_title[i]) + ylab(NULL) + xlab(NULL)
    if(seg){
      if(intg){
        tick_width <- c(round((y_max_low - miny)/2.5), round((maxy - y_min_high)/2.5))
      }else{
        tick_width <- c(floor_dec((y_max_low - miny)/2.5, 1), round((maxy - y_min_high)/2.5))
      }
      pp <- gg.gap(plot = pp, ylim = c(miny, maxy),
                   tick_width = tick_width,
                   segments = c(y_max_low, y_min_high))+
        theme(plot.margin = unit(c(1, 0, 1, 0), "lines"), legend.position = 'bottom')
    }
    pics[[i]] <- pp
  }
  if(length(n_cols) <= 4){
    ans <- ggarrange(plotlist = pics, ncol = length(n_cols), nrow = 1,
                     common.legend = TRUE, legend = 'bottom')
  }else{
    ans <- ggarrange(plotlist = pics, ncol = 4, nrow = ceiling(length(n_cols)/4),
                     common.legend = TRUE, legend ='bottom')
  }
  ans <- annotate_figure(ans, left = text_grob(y_label, rot = 90, size = 17, hjust = hjust_))
  ans
}

## summarise data information / three diversity values / evenness / completeness table
summa_info <- function(myout, myinfo, myeven) {
  myinfo <- myinfo %>% select(Year, n, SC, f1, f2) %>% as_tibble()
  myeven$goalSC <- as.character(myeven$goalSC)
  ans <- myeven %>% select(-Assemblage) %>%
    mutate(goalSC1 = gsub(goalSC, pattern = 'Coverage = ', replacement = '')) %>% 
    mutate(goalSC2 = as.numeric(gsub(goalSC1, pattern = '%', replacement = ''))/100) %>%
    mutate(even2 = round(even, 2)) %>% 
    select(Year = year, goalSC = goalSC2, even = even2) %>% 
    reshape2::dcast(data = ., formula = Year ~ goalSC, value.var = 'even') %>% as_tibble() %>% 
    left_join(x = myinfo, y = ., by = 'Year')
  #SC,q=0
  ans <- myout %>% filter(type == 'qD', type2 == 'Taxonomic', Order.q == 0) %>% 
    mutate(Year = year + 1980, goalSC1 = gsub(goalSC, pattern = 'Coverage = ', replacement = '')) %>% 
    select(Year,goalSC = goalSC1, qD = qD) %>% reshape2::dcast(data = ., formula = Year ~ goalSC, value.var = 'qD') %>% 
    as_tibble() %>% mutate(SC0 = round(Observed/Asymptotic, 4)) %>% select(Year, SC0) %>% 
    left_join(x = ans, y = ., by = 'Year') %>%
    select(Year, n, SC0, SC, f1, f2, E_min = as.character(goalSC[1]), E_max = as.character(goalSC[2]))
  #qD
  myout <- myout %>% filter(type == 'qD') %>% 
    mutate(Year = year + 1980, goalSC1 = gsub(goalSC, pattern = 'Coverage = ', replacement = '')) %>% 
    select(Year, Order.q, goalSC = goalSC1, qD = qD,type2)
  myout$goalSC <- factor(myout$goalSC, levels = c('Observed', 'Asymptotic', goalSC), labels = c('Observed', 'Asymptotic', 'C_min', 'C_max'))
  myout$type2 <- factor(myout$type2, levels = unique(myout$type2))
  ans <- reshape2::dcast(data = myout, formula = Year ~ type2 + Order.q + goalSC, value.var = 'qD') %>% as_tibble() %>% 
    left_join(x = ans, y = ., by = 'Year')
  ans[, -(1:9)] <- round(ans[, -(1:9)], 2)
  return(ans)
}

## plot coverage along with each time point.
# Arguments :
# 'seg' is the logical setting whether figure has segment or not.
cvrg_plot <- function(mycvrg, seg = FALSE) {
  cvrg_p <- ggplot(data = cvrg) +
    geom_line(aes(x = Year, y = Completeness, col = Order.q), lty = 2) +
    geom_point(aes(x = Year, y = Completeness, col = Order.q, shape = Order.q), size = 4) + theme_bw() +
    theme(text = element_text(size = 17), legend.text = element_text(size = 17), legend.position = 'bottom',
          legend.key.size = unit(3, "line"), legend.margin = margin(-0.3, 0.2, -0.3, 0.2, "cm")) +
    scale_color_manual(values = c('#1F78B4', '#E7298A')) +
    scale_shape_manual(values = c(17, 16))
  if (seg) {
    floor_dec <- function(x, dig = 1) round(x - 5*10^(-dig-1), dig)
    ceiling_dec <- function(x, dig = 1) round(x + 5*10^(-dig-1), dig)
    maxy <- 1
    miny <- floor_dec(min(mycvrg$Completeness), 2)
    y_max_low <- ceiling_dec(max(mycvrg$Completeness[mycvrg$Order.q == 0], na.rm = T), 2)
    y_min_high <- floor_dec(min(mycvrg$Completeness[mycvrg$Order.q == 1], na.rm = T), 3)
    tick_width <- c(floor_dec((y_max_low - miny)/2.5, 2), round((maxy - y_min_high)/2.5, 3))
    
    cvrg_p <- gg.gap(plot = cvrg_p, ylim = c(miny, maxy), tick_width = tick_width,
                     segments = c(y_max_low, y_min_high))
  }
  return(cvrg_p)
}


######## abundance data analysis -------------------------------------------------------------------
q = c(0, 1, 2)

## read.data
dataset <- read_csv("Fish monthly data.csv")

## window width (time width in each group)
wd = 1
mydata <- abun(dataset, wd)


## overall abundance / rare group
method = 'overall'    # 'overall' or 'rare'
addabu <- FALSE       # 'FALSE' means only using rare species withou abundant species

## Calculate basic Data information
if (method == 'overall') {
  infos <- iNEXT3D:::DataInfo(x = rbind(mydata$abun, mydata$rare), datatype = 'abundance') %>% 
    mutate(Year = as.numeric(substr(Assemblage, 3, 8))) %>% select(-Assemblage) 
} else {
  infos <- iNEXT3D:::DataInfo(x = mydata$rare, datatype = 'abundance') %>% 
    mutate(Year = as.numeric(substr(Assemblage, 3, 8))) %>% select(-Assemblage) 
}

## summary of sample coverage (q = 1)
summary(infos$SC)

## Calculate Cminimum and Cmaximum
if (method == 'overall') {
  cmax <- apply(rbind(mydata$abun, mydata$rare), 2, function(x) iNEXT3D:::Chat.Ind(x, 2*sum(x))) %>% 
    min %>% round(., 4)
  cmin <- apply(rbind(mydata$abun, mydata$rare), 2, function(x) iNEXT3D:::Chat.Ind(x, sum(x))) %>% 
    min %>% round(., 4)
} else {
  cmax <- apply(mydata$rare, 2, function(x) iNEXT3D:::Chat.Ind(x, 2*sum(x))) %>% 
    min %>% round(., 4)
  cmin <- apply(mydata$rare, 2, function(x) iNEXT3D:::Chat.Ind(x, sum(x))) %>% 
    min %>% round(., 4)
}
goalSC <- c(cmin, cmax)

######## Taxonomic -------------------------------------------------------------------
## Calculate Taxonomic diversity under Cmin and Cmax
if (method == 'overall') {
  out_td <- estimate3D(data = rbind(mydata$abun,mydata$rare), class = 'TD', q = q, datatype = 'abundance', 
                       base = 'coverage', level = goalSC, nboot = 0) 
} else {
  out_td <- estimate3D(data = mydata$rare, class = 'TD', q = q, datatype = 'abundance', 
                       base = 'coverage', level = goalSC, nboot = 0) 
}

out_td <- out_td %>% select(Assemblage, Order.q, qD, goalSC)
out_td$goalSC <- as.character(out_td$goalSC)

 
## Calculate asymptotic Taxonomic diversity and empirical Taxonomic diversity 
if (method == 'overall') {
  out_td_emp <- Obs3D(data = rbind(mydata$abun, mydata$rare), class = 'TD', q = q, datatype = 'abundance', nboot = 0) %>% 
    select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed') 
  out_td_est <- Asy3D(data = rbind(mydata$abun,mydata$rare), q = q, datatype = 'abundance', nboot = 0) %>% 
    select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic') 
} else {
  out_td_emp_abun <- Obs3D(data = mydata$abun, class = 'TD', q = q, datatype = 'abundance', nboot = 0) %>% 
    select(Assemblage,Order.q,qD) 
  out_td_emp_abun$Assemblage <- as.character(out_td_emp_abun$Assemblage)
  if (addabu == FALSE) {
    out_td_emp_abun$qD <- 0
  }
  out_td <- left_join(x = out_td, y = out_td_emp_abun, by = c("Assemblage", 'Order.q')) %>% mutate(qD = qD.x + qD.y) %>% 
    select(-qD.x, -qD.y) %>% select(Assemblage, Order.q, qD, goalSC)
  
  out_td_emp_rare <- Obs3D(data = mydata$rare, class = 'TD', q = q, datatype = 'abundance', nboot = 0) %>% 
    select(Assemblage, Order.q = Order.q, qD) %>% mutate(goalSC = 'Observed') 
  out_td_emp <- out_td_emp_rare
  out_td_emp$qD <- out_td_emp$qD + out_td_emp_abun$qD
  
  out_td_est_rare <- Asy3D(data = mydata$rare, class = 'TD', q = q, datatype = 'abundance', nboot = 0) %>% 
    select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic')
  out_td_est <- out_td_est_rare
  out_td_est$qD <- out_td_est$qD + out_td_emp_abun$qD
  
  rm(out_td_emp_abun); rm(out_td_emp_rare); rm(out_td_est_rare)
}


## Plot Evenness
even <- out_td %>% filter(Order.q == 1 | Order.q == 0) %>%
  filter(goalSC == as.character(goalSC[1]) | goalSC == as.character(goalSC[2])) %>% 
  group_by(Assemblage, goalSC) %>% 
  summarise(even = (qD[Order.q == 1] - 1)/(qD[Order.q==0] - 1)) %>% ungroup() %>% 
  mutate(year = as.numeric(substr(Assemblage, 3, nchar(Assemblage))))
even$goalSC[even$goalSC == as.character(goalSC[1])] <- paste0('Coverage = ', goalSC[1]*100, '%')
even$goalSC[even$goalSC == as.character(goalSC[2])] <- paste0('Coverage = ', goalSC[2]*100, '%')
even$goalSC <- factor(even$goalSC,levels = unique(even$goalSC)) 
evenp <- ggplot(data = even) + geom_line(aes(x = year, y = even, color = goalSC), lty = 2) +
  geom_point(aes(x = year, y = even, col = goalSC, shape = goalSC), size = 4) +
  theme_bw() + xlab('Year') +
  theme(legend.position = 'bottom', legend.title = element_blank(),text = element_text(size = 17),
        legend.key.size = unit(3, "line"), legend.margin = margin(-0.3, 0.2, -0.3, 0.2, "cm")) +
  scale_color_manual(values = c('#1F78B4', '#E7298A')) + ylab(label = 'Evenness')
evenp


## Check significance
out_td <- rbind(out_td, out_td_emp, out_td_est)
out_td <- find_ord(myout = out_td)

## plot q-profile Taxonomic diversity under each coverage.
pic_td <- pic_qs(myout = out_td, n_cols = 1:(length(goalSC)+2), y_label = 'Taxonomic diversity')
pic_td

######## Phylogenetic -------------------------------------------------------------------
## read data and phylo tree
phylotr <- read.tree("Fish PhyTree.txt")
phylotr$tip.label[which(phylotr$tip.label == "Gadus_morhua_-species_in_domain_Eukaryota")] <- 'Gadus_morhua'
phylotr$tip.label[which(phylotr$tip.label == "Trigloporus_lastoviza")] <- 'Chelidonichthys_lastoviza'
sum(phylotr$tip.label %in% rownames(rbind(mydata$abun, mydata$rare)))
phyclust::get.rooted.tree.height(phylotr)    ## get the tree height

## Calculate Phylogenetic diversity under Cmin and Cmax
if (method == 'overall') {
  out_pd <- estimate3D(data = rbind(mydata$abun,mydata$rare), class = 'PD', tree = phylotr, datatype = 'abundance', q = q, PDtype = 'meanPD',
                       level = goalSC, nboot = 0, reftime = 1) %>% select(Assemblage, Order.q, qPD, goalSC)
} else {
  out_pd <- estimate3D(data = mydata$rare, tree = phylotr, class = 'PD', datatype = 'abundance', q = q, PDtype = 'meanPD',
                       level = goalSC, nboot = 0, reftime = 1) %>% select(Assemblage, Order.q, qPD, goalSC)
}
out_pd$goalSC <- as.character(out_pd$goalSC)

## Calculate asymptotic Phylogenetic diversity and empirical Phylogenetic diversity 
if (method == 'overall') {
  out_pd_emp <- Obs3D(data = rbind(mydata$rare, mydata$abun), class = 'PD', tree = phylotr, datatype = 'abundance', PDtype = 'meanPD',
                      q = q, nboot = 0, reftime = 1)  %>%
    select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Observed') 
  out_pd_est <- Asy3D(data = rbind(mydata$abun,mydata$rare),class = 'PD', tree = phylotr, datatype = 'abundance', PDtype = 'meanPD',
                      q = q, nboot = 0, reftime = 1) %>% 
    select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic') 
} else {
  out_pd_emp_abun <- Obs3D(data = mydata$abun, class = 'PD', tree = phylotr, datatype = 'abundance', PDtype = 'meanPD',
                           q = q, nboot = 0, reftime = 1) %>%
    select(Assemblage, Order.q, qPD) 
  if (addabu == FALSE) {
    out_pd_emp_abun$qPD <- 0
  }
  out_pd <- left_join(x = out_pd,y = out_pd_emp_abun, by = c("Assemblage",'Order.q')) %>% mutate(qPD = qPD.x + qPD.y) %>% 
    select(-qPD.x, -qPD.y) %>% select(Assemblage, Order.q, qPD, goalSC)
  
  out_pd_emp_rare <- Obs3D(data = mydata$rare, class = 'PD', tree = phylotr, datatype = 'abundance',
                           PDtype = 'meanPD', q = q, nboot = 0, reftime = 1) %>%
    select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Observed') 
  out_pd_emp <- out_pd_emp_rare
  out_pd_emp$qPD <- out_pd_emp$qPD + out_pd_emp_abun$qPD
  
  out_pd_est_rare <- Asy3D(data = mydata$rare, class = 'PD', tree = phylotr, datatype = 'abundance', 
                           PDtype = 'meanPD', q = q, nboot = 0, reftime = 1) %>% 
    select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic') 
  out_pd_est <- out_pd_est_rare
  out_pd_est$qPD <- out_pd_est$qPD + out_pd_emp_abun$qPD
  
  rm(out_pd_emp_abun); rm(out_pd_emp_rare)
}


## Check significance
out_pd <- rbind(out_pd, out_pd_emp, out_pd_est)
out_pd <- find_ord(myout = out_pd)

## plot q-profile Phylogenetic diversity under each coverage.
pic_pd <- pic_qs(myout = out_pd, n_cols = 1:(length(goalSC)+2), y_label = 'Phylogenetic diversity', nbreak = 6)
pic_pd

######## Functional -------------------------------------------------------------------
## specify the number of traits for distance matrix: 5 or 12
traits_used <- 12   ## User can select 5 or 12

## read traits
traits <- read_csv("Fish Traits.csv")
if (traits_used == 5) {
  traits <- traits[, c('Species', 'meanWeight hinkley', 'Lm', 'Troph', 'TempPrefMean', 'PositionWaterColumn')]
} else if (traits_used == 12) {
  traits <- traits[, c(1, 2, 3, 7:16)]
}

## Turn traits into pairwise distance matrix in gower distance between zero and 1.
traits <- as.data.frame(traits)
rownames(traits) <- traits$Species; traits <- traits[, -1]
traits <- traits[match(rownames(rbind(mydata$abun,mydata$rare)), rownames(traits)), ]

## Check the order of species in traits table and data
sum(rownames(traits) == rownames(rbind(mydata$abun, mydata$rare)))

## Turn character traits into factor and generate distance matrix
for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[,i] <- factor(traits[,i], levels = unique(traits[,i]))
}
disM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


## Calculate Functional diversity under Cmin and Cmax
if (method == 'overall') {
  out_fd <- estimate3D(data = rbind(mydata$abun, mydata$rare), class = 'AUC', distM = disM, q = q, 
                       level = goalSC, nboot = 0, datatype = 'abundance') %>% 
    select(Assemblage, Order.q, qFD = qAUC, goalSC)
} else {
  rare_index <- match(rownames(mydata$rare), rownames(disM))
  disM_rare <- disM[rare_index, rare_index]
  out_fd <- estimate3D(data = mydata$rare, class = 'AUC', distM = disM_rare, q = q,
                       level = goalSC, nboot = 0, datatype = 'abundance') %>% 
    select(Assemblage, Order.q, qFD = qAUC, goalSC)
  rm(rare_index); rm(disM_rare)
}
out_fd$goalSC <- as.character(out_fd$goalSC)

## Calculate asymptotic Functional diversity and empirical Functional diversity 
if (method == 'overall') {
  out_fd_emp <- Obs3D(data = rbind(mydata$abun,mydata$rare), class = 'AUC',
                      distM = disM, q = q, datatype = 'abundance', nboot = 0) %>%
    select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed') 
  
  out_fd_est <- Asy3D(data = rbind(mydata$abun, mydata$rare), class = 'AUC',
                      distM = disM, datatype = 'abundance', q, nboot = 0) %>% 
    select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic')
} else {
  abun_index <- match(rownames(mydata$abun), rownames(disM))
  rare_index <- match(rownames(mydata$rare), rownames(disM))
  disM_abun <- disM[abun_index, abun_index]
  disM_rare <- disM[rare_index, rare_index]
  out_fd_emp_abun <- Obs3D(data = mydata$abun, class = 'AUC', distM = disM_abun, q = q,
                           datatype = 'abundance', nboot = 0) %>% ungroup() %>% 
    select(Assemblage, Order.q, qFD = qAUC)
  
  if (addabu == FALSE) out_fd_emp_abun$qFD <- 0
  
  out_fd <- left_join(x = out_fd,y = out_fd_emp_abun, by = c("Assemblage", 'Order.q')) %>% 
    mutate(qFD = qFD.x + qFD.y) %>% select(-qFD.x, -qFD.y) %>% select(Assemblage, Order.q, qFD, goalSC)
  out_fd_emp_rare <- Obs3D(data = mydata$rare, class = 'AUC',distM = disM_rare,q = q,
                           datatype = 'abundance', nboot = 0) %>% ungroup() %>% 
    select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed') 
  out_fd_emp <- out_fd_emp_rare
  out_fd_emp$qFD <- out_fd_emp$qFD + out_fd_emp_abun$qFD
  out_fd_est_rare <- Asy3D(data = mydata$rare, class = 'AUC',distM = disM_rare, q = q,
                           datatype = 'abundance', nboot = 0) %>% ungroup() %>% 
    select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic') %>% ungroup()
  out_fd_est <- out_fd_est_rare
  out_fd_est$qFD <- out_fd_est$qFD + out_fd_emp_abun$qFD
  
  rm(abun_index); rm(rare_index)
  rm(disM_abun); rm(disM_rare)
  rm(out_fd_emp_abun); rm(out_fd_emp_rare)
  rm(out_fd_est_rare)
}

## Check significance
out_fd <- rbind(out_fd, out_fd_emp, out_fd_est)
out_fd <- find_ord(myout = out_fd)

## plot q-profile Functional diversity under each coverage.
pic_fd <- pic_qs(myout = out_fd,n_cols = 1:(length(goalSC)+2), y_label = 'Functional diversity', nbreak = 7)
pic_fd

######## summarise three diversity as the appendix table -------------------------------------------------------------------
out_td <- out_td %>% mutate(type2 = 'Taxonomic')
out_pd <- out_pd %>% mutate(type2 = 'Phylogenetic')
out_fd <- out_fd %>% mutate(type2 = 'Functional')
out_tpfd <- rbind(out_td, out_pd, out_fd)
infos <- summa_info(myout = out_tpfd, myinfo = infos, myeven = even)
# write.csv(infos, file = 'infos.csv', row.names = F)

y_label <- 'Diversity'
if (method == 'overall') {
  pic_tpfd0 <- pic_tpf(myout = out_tpfd, qfilter = 0, n_cols = 1:4, y_label, seg = T, intg = T, max_2nd = F)
} else if (method == 'sep') {
  pic_tpfd0 <- pic_tpf(myout = out_tpfd, qfilter = 0, n_cols = 1:4, y_label, seg = F, intg = T, max_2nd = F)
}
pic_tpfd0

pic_tpfd1 <- pic_tpf(myout = out_tpfd, qfilter = 1, n_cols = 1:4, y_label, seg = F, intg = F)
pic_tpfd1

pic_tpfd2 <- pic_tpf(myout = out_tpfd, qfilter = 2, n_cols = 1:4, y_label, seg = F, intg = F)
pic_tpfd2

## plot coverage
cvrg <- infos[, c('Year', 'SC', 'SC0')] %>% melt(id.vars = c('Year'), value.name = 'Completeness') %>% 
  mutate(Order.q = factor(ifelse(variable == 'SC', 1, 0))) %>% select(-variable) %>% as_tibble()

if (method == 'overall') {
  cvrg_p <- cvrg_plot(cvrg, seg = TRUE)
} else if (method == 'sep') {
  cvrg_p <- cvrg_plot(cvrg, seg = FALSE)
}
cvrg_p

######## incidence data analysis -------------------------------------------------------------------
q = c(0, 1, 2)

## read.data
dataset <- read_csv("Fish monthly data.csv")

## window width (time width in each group)
wd = 3
mydata <- inci_freq(dataset, wd)
nT <- mydata[1,]

## overall abundance / rare group
method = 'overall'    # 'overall' or 'rare'

## Calculate basic Data information
infos <- iNEXT3D:::DataInfo(x = mydata, datatype = 'abundance') %>% 
  mutate(Year = as.numeric(substr(Assemblage, 3, 8))) %>% select(-Assemblage) 

## summary of sample coverage (q = 1)
summary(infos$SC)

## Calculate Cminimum and Cmaximum
cmax <- apply(mydata, 2, function(x) iNEXT3D:::Chat.Sam(x, 2*x[1])) %>% min %>% round(., 4)
cmin <- apply(mydata, 2, function(x) iNEXT3D:::Chat.Sam(x, x[1])) %>% min %>% round(., 4)
goalSC <- c(cmin,cmax)

######## Taxonomic -------------------------------------------------------------------
## Calculate Taxonomic diversity under Cmin and Cmax
out_td <- estimate3D(data = mydata, class = 'TD', q = q, datatype = 'incidence_freq', base = 'coverage', 
                     level = goalSC, nboot = 0) 

out_td <- out_td %>% select(Assemblage, Order.q, qD, goalSC)
out_td$goalSC <- as.character(out_td$goalSC)

## Calculate asymptotic Taxonomic diversity and empirical Taxonomic diversity 
out_td_emp <- Obs3D(data = mydata, class = 'TD', q = q, datatype = 'incidence_freq', nboot = 0) %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Observed') 
out_td_est <- Asy3D(data = mydata, class = 'TD', q = q, datatype = 'incidence_freq', nboot = 0) %>% select(Assemblage, Order.q, qD) %>% mutate(goalSC = 'Asymptotic') 

## Plot Evenness
even <- out_td %>% filter(Order.q == 1 | Order.q == 0 ) %>%
  filter(goalSC == as.character(goalSC[1])|goalSC == as.character(goalSC[2])) %>% 
  group_by(Assemblage, goalSC) %>% 
  summarise(even = (qD[Order.q == 1] - 1)/(qD[Order.q == 0] - 1)) %>% ungroup() %>% 
  mutate(year = as.numeric(substr(Assemblage, 3, nchar(Assemblage))))
even$goalSC[even$goalSC == as.character(goalSC[1])] <- paste0('Coverage = ', goalSC[1]*100, '%')
even$goalSC[even$goalSC == as.character(goalSC[2])] <- paste0('Coverage = ', goalSC[2]*100, '%')
even$goalSC <- factor(even$goalSC,levels = unique(even$goalSC)) 
evenp <- ggplot(data = even) + geom_line(aes(x = year, y = even, color = goalSC), lty = 2)+
  geom_point(aes(x = year, y = even,col = goalSC, shape = goalSC), size = 4)+
  theme_bw() + xlab('Year')+
  theme(legend.position = 'bottom', legend.title = element_blank(), text = element_text(size = 17),
        legend.key.size = unit(3, "line"),legend.margin = margin(-0.3, 0.2, -0.3, 0.2, "cm"))+
  scale_color_manual(values = c('#1F78B4', '#E7298A')) + ylab(label = 'Evenness')
evenp


## Check significance
out_td <- rbind(out_td, out_td_emp, out_td_est)
out_td <- find_ord(myout = out_td)

## plot q-profile Taxonomic diversity under each coverage.
pic_td <- pic_qs(myout = out_td,n_cols = 1:(length(goalSC) + 2), y_label = 'Taxonomic diversity')
pic_td

######## Phylogenetic -------------------------------------------------------------------
## read data and phylo tree
phylotr <- read.tree("Fish PhyTree.txt")
phylotr$tip.label[which(phylotr$tip.label == "Gadus_morhua_-species_in_domain_Eukaryota")] <- 'Gadus_morhua'
phylotr$tip.label[which(phylotr$tip.label == "Trigloporus_lastoviza")] <- 'Chelidonichthys_lastoviza'

## Calculate Phylogenetic diversity under Cmin and Cmax
out_pd <- estimate3D(data = inci_raw(dataset, wd), class = 'PD', tree = phylotr, datatype = 'incidence_raw', q = q, PDtype = 'meanPD', nT = nT,
                     level = goalSC, nboot = 0, reftime = 1) %>% select(Assemblage, Order.q, qPD, goalSC)
out_pd$goalSC <- as.character(out_pd$goalSC)

## Calculate asymptotic Phylogenetic diversity and empirical Phylogenetic diversity 
out_pd_emp <- Obs3D(data = inci_raw(dataset,wd), class = 'PD', tree = phylotr, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                    q = q, nboot = 0, reftime = 1) %>%
  dplyr::select(Assemblage, Order.q, qPD ) %>% mutate(goalSC = 'Observed')

out_pd_est <- Asy3D(data = inci_raw(dataset, wd), class = 'PD', tree = phylotr, datatype = 'incidence_raw', PDtype = 'meanPD', nT = nT,
                    q = q, nboot = 0, reftime = 1) %>% 
  dplyr::select(Assemblage, Order.q, qPD) %>% mutate(goalSC = 'Asymptotic')


## Check significance
out_pd <- rbind(out_pd, out_pd_emp, out_pd_est)
out_pd <- find_ord(myout = out_pd)

## plot q-profile Phylogenetic diversity under each coverage.
pic_pd <- pic_qs(myout = out_pd,n_cols = 1:(length(goalSC)+2), y_label = 'Phylogenetic diversity', nbreak = 6)
pic_pd

######## Functional -------------------------------------------------------------------
## specify the number of traits for distance matrix: 5 or 12
traits_used <- 12   ## User can select 5 or 12

## read traits
traits <- read_csv("Fish Traits.csv")
if (traits_used == 5) {
  traits <- traits[, c('Species', 'meanWeight hinkley', 'Lm', 'Troph', 'TempPrefMean', 'PositionWaterColumn')]
} else if (traits_used == 12) {
  traits <- traits[, c(1, 2, 3, 7:16)]
}

## Turn traits into pairwise distance matrix in gower distance between zero and 1.
traits <- as.data.frame(traits)
rownames(traits) <- traits$Species; traits <- traits[, -1]
traits <- traits[match(rownames(mydata), rownames(traits)), ]

## Check the order of species in traits table and data
sum(rownames(traits) == rownames(mydata))

## Turn character traits into factor and generate distance matrix
for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[, i] <- factor(traits[,i], levels = unique(traits[, i]))
}
disM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


## Calculate Functional diversity under Cmin and Cmax
out_fd <- estimate3D(data = mydata, class = 'AUC', distM = disM[-1,-1], q = q,
                     level = goalSC, nboot = 0, datatype = 'incidence_freq') %>% 
  select(Assemblage, Order.q, qFD = qAUC, goalSC)

out_fd$goalSC <- as.character(out_fd$goalSC)


## Calculate asymptotic Functional diversity and empirical Functional diversity 
out_fd_emp <- Obs3D(data = mydata, class = 'AUC', distM = disM[-1, -1], datatype = 'incidence_freq', 
                    q = q, nboot = 0) %>% 
  select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Observed')

out_fd_est <- Asy3D(data = mydata, class = 'AUC', distM = disM[-1,-1], datatype = 'incidence_freq',
                    q = q, nboot = 0) %>% 
  select(Assemblage, Order.q, qFD = qAUC) %>% mutate(goalSC = 'Asymptotic')


## Check significance
out_fd <- rbind(out_fd, out_fd_emp, out_fd_est)
out_fd <- find_ord(myout = out_fd)

## plot q-profile Functional diversity under each coverage.
pic_fd <- pic_qs(myout = out_fd,n_cols = 1:(length(goalSC)+2), y_label = 'Functional diversity', nbreak = 7)
pic_fd

######## summarise TPFD as the appendix table -------------------------------------------------------------------
out_td <- out_td %>% mutate(type2 = 'Taxonomic')
out_pd <- out_pd %>% mutate(type2 = 'Phylogenetic')
out_fd <- out_fd %>% mutate(type2 = 'Functional')
out_tpfd <- rbind(out_td, out_pd, out_fd)
infos <- summa_info(myout = out_tpfd,myinfo = infos, myeven = even)
# write.csv(infos, file = 'infos.csv', row.names = F)

y_label <- 'Diversity'
pic_tpfd0 <- pic_tpf(myout = out_tpfd, qfilter = 0, n_cols = 1:4, y_label, seg = T, intg = T, max_2nd = F)
pic_tpfd0

pic_tpfd1 <- pic_tpf(myout = out_tpfd, qfilter = 1, n_cols = 1:4, y_label, seg = F, intg = F)
pic_tpfd1

pic_tpfd2 <- pic_tpf(myout = out_tpfd,qfilter = 2,n_cols = 1:4, y_label, seg = F, intg = F)
pic_tpfd2

## plot coverage
cvrg <- infos[,c('Year', 'SC', 'SC0')] %>% melt(id.vars = c('Year'), value.name = 'Completeness') %>% 
  mutate(Order.q = factor(ifelse(variable == 'SC', 1, 0))) %>% select(-variable) %>% as_tibble()

if (method == 'overall') {
  cvrg_p <- cvrg_plot(cvrg, seg = TRUE)
} else if (method == 'sep') {
  cvrg_p <- cvrg_plot(cvrg, seg = FALSE)
}

cvrg_p



