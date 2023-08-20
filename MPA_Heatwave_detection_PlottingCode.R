# April 2023
# Author: Jess K Hopf
# to be run with data generated in MPA_Heatwave_detection_v7.m

# Load packages
library(tidyverse); library(viridisLite); library(scales); library(Ecfun)
library(ggridges)
 
# Clear environment
rm(list=ls())

# load required data
# this is the time is takes to reach high detectability in the absence of 
# distrubance for the baseline and the sensitivity analysis scenarios
# df name is T2AUCnoDistData
load("Times2AUCNoDist_20230508.Rda")


# Functions ---------------

# clean up dataframe
  clean_func <- function(Results_in){
        
      Results <- pivot_longer(Results_in, cols = IO:BACI, 
                              names_to = "Sample_method",
                              values_to = "AUC")
      
      # change Dist_Impact to proportion removed with disturbance
      # (keeps the sign for sign pop scenarios, sign sign of 0 = pos)
      Results <- Results %>%
        mutate(Prop_gone = (1-abs(Dist_Impact))*Ecfun::sign(Dist_Impact, zero=1)) 
      
      # clean table
      Results$Prop_gone <- factor(Results$Prop_gone) # factor(round(Results$Prop_gone, 2))
      Results$Sample_method <- ordered(Results$Sample_method, 
                                          levels = c("IO", "BA", "BACI"))
      
      
      # if AUC results are <0.5 do NA
      Results_05 <- Results %>%
                    mutate(AUC = ifelse(AUC < 0.5, NA, AUC))
      
      return(Results_05)
  }


# time taken to reach AUC >=0.8 (high detectability)
  time2HD_func <- function(Results_05, scen){
    # less values of prop_gone
    PG_vec <- seq(0,1,0.1)  
    
    # add time to HD without disturbance
    # (it takes a certain amount of time to reach high detectability in the absence of disturbance)
    # calculate if the disturbance is before or after it reaches HD, by sampling method 
    Results <- left_join(Results_05,
                         T2AUCnoDistData %>% filter(scenSA == scen) %>% select(-scenSA), 
                         by = "Sample_method") %>% 
                mutate(isbeforeT2HD = Yr_Dist_Start < T2HDwoDist,
                       isSampb4T2HD = Yr_Sampled < T2HDwoDist)
    
    #     1. Disturbance happens before it reaches high detect
    #         - this is the first time the AUC is reached (might extend AUC out)
    Time2AUC1 <- Results %>% filter(Prop_gone %in% PG_vec, 
                                    isbeforeT2HD == TRUE) %>% 
      arrange(Yr_Dist_Start, Dist_Length, Sample_method, Prop_gone, Yr_Sampled) %>% 
      group_by(Yr_Dist_Start, Dist_Length, Sample_method, Prop_gone) %>% 
      filter(AUC >= 0.80) %>% slice_head(n=1) %>% 
      mutate(Prop_gone = ordered(Prop_gone, levels = PG_vec),
             BAdist = 'Before') 
    
    #     2. Distrubance happens after high detect is reached 
    #         - in this case its the second time AUC is reached (if affected at all)
    #         - will take out first 4 years & all years before disturbance ends
    Time2AUC2 <- Results %>% filter(Prop_gone %in% PG_vec, 
                                    isSampb4T2HD == FALSE,
                                    isbeforeT2HD == FALSE,
                                    Samp_time >0) %>% 
      arrange(Yr_Dist_Start, Dist_Length, Sample_method, Prop_gone, Yr_Sampled) %>% 
      group_by(Yr_Dist_Start, Dist_Length, Sample_method, Prop_gone) %>% 
      filter(AUC >= 0.80) %>% slice_head(n=1) %>% 
      mutate(Prop_gone = ordered(Prop_gone, levels = PG_vec),
             BAdist = 'After')
    
    # bind
    # last command creates a col to assign text colour in heatmaps
    Time2AUC <- bind_rows(Time2AUC1, Time2AUC2) %>% 
      mutate(Prop_goneAlpha = ordered(as.numeric(as.character(Prop_gone)),levels = PG_vec),
             Dist_Length = ordered(Dist_Length, levels = c(1,3,5)),
             Yr_SampledMax = max(Yr_Sampled, Yr_Dist_End),
             Samp_timeMax = max(0, Samp_time),
             HMtextcol = Samp_timeMax>9) 
    
    return(Time2AUC)
  }  
  
  
  
# FIGURES ------------  
  
# Baseline main figures
  # read in file
  # choose either adults or juvs/ open or closed
  Results_in <- read_csv("Outputs/2023-03-09_V7results_open_adults_BASELINE.csv")
  # Results_in <- read_csv("Outputs/2023-03-09_V7results_closed_adults_BASELINE.csv")
  # Results_in <- read_csv("Outputs/2023-03-06_V7results_open_juvs_BASELINE.csv")
  Results_05 <- clean_func(Results_in)

  
# Contour lines (Main fig)
  ggplot(Results_05 %>% filter(Yr_Dist_Start %in% c(0,4,8))) +
    facet_grid(Dist_Length~Yr_Dist_Start , labeller=label_both) +
    geom_contour(aes(x = Yr_Sampled, y = as.numeric(as.character(Prop_gone)),
                     z = AUC, color = stat(level), size = Sample_method), 
                 breaks = c(0.6,0.8)) +  # breaks = c(0.6,0.8)) +
    scale_colour_viridis_c(direction = -1, option = "C", limits = c(0.5,1)) +
    scale_size_manual(values = c(0.5,2,4)) +
    # scale_linetype_manual(values = c("dotted", "dotdash", "solid")) +
    theme_light() +
    scale_x_continuous(breaks = seq(0,25,by=4), 
                       minor_breaks = seq(0,25,by=2), 
                       name = 'Time since MPA est.') +
    scale_y_continuous(name = 'Proportion reduction in abundance (0 = no dist.)') 
  
  
  # Save (NOTE this can't handle transparencies)
  # ggsave(filename = "Contour_Adults_Both_BASELINE.eps",
  #        width = 19, height = 12, units = "cm",
  #        device = "eps")
  

# time to high detectability 
Time2AUCBL <- time2HD_func(Results_05, "BL")
  
  # some stats
  Time2AUCBL %>% filter(Yr_Dist_Start < 4) %>% ungroup() %>% 
    summarise(mean = mean(Samp_timeMax),
              median = median(Samp_timeMax))
  
  Time2AUCBL %>% filter(Yr_Dist_Start > 3) %>% ungroup() %>% 
    summarise(mean = mean(Samp_timeMax),
              median = median(Samp_timeMax))
  

  # scatter plot (Main Fig)
  ggplot(Time2AUCBL %>% filter(Yr_Dist_Start %in% c(0,4,8)),
         aes(y = Samp_timeMax, x = Dist_Length, color = Prop_gone)) +
    geom_point(alpha = 0.7, size = 2.5, position = position_dodge(0.1)) +
    facet_grid(Sample_method~Yr_Dist_Start, labeller = label_both) +
    scale_color_viridis_d(direction = 1, option = "C", end = 0.8) +
    scale_y_continuous(breaks = pretty_breaks(),
                       name = 'Time to High (0.8) Detectability after disturbance ends') +
    scale_x_discrete(name = "Length of Disturbance (Yrs)") +
    theme_light()
  
  
  # heatmap (appendix)
  # single dataset
  ggplot(Time2AUCBL %>% filter(Yr_Dist_Start %in% c(0,4,8)) %>% 
           mutate("MPA age" = Yr_Dist_Start),
           aes(y = Prop_gone,
               x = Dist_Length)) +
      geom_tile(aes(fill = factor(Samp_timeMax, levels = c(0:14)))) +
      geom_text(aes(label = Samp_timeMax, colour = HMtextcol)) + 
      facet_grid(Sample_method~`MPA age`, 
                 labeller = labeller(.cols = label_both)) +
      scale_fill_viridis_d(option = "E", drop = FALSE,
                           na.value = "#FDE725FF",
                           labels = c(as.character(0:14), "15+"),
                           name = "Time to high detectability (yrs)",
                           guide = guide_legend(title.position = "right",
                             title.theme = element_text(angle = -90, 
                                                        hjust = 0.5,
                                                        vjust = 2))) +
      scale_color_manual(values = c("white", "black"), guide="none") +
      scale_y_discrete(name = 'Proportion reduction in abundance') +
      scale_x_discrete(name = "Length of Disturbance (Yrs)") +
      theme_light() 

  
  # Save (NOTE this can't handle transparencies)
  # ggsave(filename = "Open_Adult_Baseline_scatter.eps",
  #        width = 20, height = 10, units = "cm",
  #        device = "eps")
  # 
  

# Other figures -----------
  
  # high detect. heatmap
  # multiple datasets (sensitivity analysis)
  
  # read in files, clean, and calculate
  # baseline
  Results_in <- read_csv("Outputs/2023-03-09_V7results_open_adults_BASELINE.csv")
  Results_05 <- clean_func(Results_in)
  Time2AUCBL <- time2HD_func(Results_05) %>% filter(Yr_Dist_Start %in% c(0,4,8))
  
  # sensitivity data (need to chose which file to load)
  Results_in <- read_csv("Outputs/2023-05-07_V7results_closed_adults_BASELINE.csv")
  Results_05 <- clean_func(Results_in)
  Time2AUCSA <- time2HD_func(Results_05) %>% filter(Yr_Dist_Start %in% c(0,4,8))
  
  
  # join data
  Time2AUCjoin <- full_join(
    Time2AUCBL %>% 
      select(Prop_gone, Dist_Length, Samp_timeMax, Sample_method,
             Yr_Dist_Start, HMtextcol),
    Time2AUCSA %>% 
      select(Prop_gone, Dist_Length, Samp_timeMax, Sample_method,
             Yr_Dist_Start, HMtextcol), 
    by = c("Prop_gone", "Dist_Length", "Sample_method", "Yr_Dist_Start"),
    suffix = c("BL","SA")) %>% 
    mutate("MPA age" = Yr_Dist_Start)
  
  # calc diff between max times
  Time2AUCjoin <- Time2AUCjoin %>%
    mutate(DeltaMaxt = Samp_timeMaxSA - Samp_timeMaxBL,
           text = paste(Samp_timeMaxSA, "(", DeltaMaxt, ")", sep = ""))
  
  
  ggplot(Time2AUCjoin,
         aes(y = Prop_gone,
             x = Dist_Length)) +
    geom_tile(aes(fill = factor(Samp_timeMaxSA, levels = c(0:14)))) +
    geom_text(aes(label = text, colour = HMtextcolSA)) +
    facet_grid(Sample_method~`MPA age`, 
               labeller = labeller(.cols = label_both)) +
    scale_fill_viridis_d(option = "E", drop = FALSE,
                         na.value = "#FDE725FF",
                         labels = c(as.character(0:14), "15+"),
                         name = "Time to high detectability (yrs)",
                         guide = guide_legend(title.position = "right",
                                              title.theme = element_text(angle = -90, 
                                                                         hjust = 0.5,
                                                                         vjust = 2))) +
    scale_color_manual(values = c("white", "black"), guide="none") +
    scale_y_discrete(name = 'Proportion reduction in abundance') +
    scale_x_discrete(name = "Length of Disturbance (Yrs)") +
    theme_light() 
  
  
  # # Save (NOTE this can't handle transparencies)
  # ggsave(filename = "JuvsBoth_AUC0.8_20by200.eps",
  #        width = 20, height = 10, units = "cm",
  #        device = "eps")
  
  

  
# Sensitivity analysis summary ----
  
  # create function to read in and clean
  SA_readclean <- function(dataname, prefix){
    Results_in <- read_csv(dataname)
    Results_05 <- clean_func(Results_in)
    Time2AUC <- time2HD_func(Results_05, prefix) %>% 
                filter(Yr_Dist_Start %in% c(0,4,8)) %>% 
                select(Prop_gone, Dist_Length, Samp_timeMax, 
                       Sample_method, Yr_Dist_Start, T2HDwoDist, BAdist) %>% 
                mutate(scen = prefix)
    
    return(Time2AUC)              }
  
  
  # read-in basline
  dfBL <- SA_readclean("Outputs/2023-03-09_V7results_open_adults_BASELINE.csv", "BL")
  
  # read-in and join all SA
  dfSAs <- rbind(SA_readclean("Outputs/2023-03-09_V7results_open_adults_BASELINE.csv", "BL"),
                 SA_readclean("Outputs/2023-03-06_V7results_open_juvs_BASELINE.csv", "BLjuvs"),
                 SA_readclean("Outputs/2023-05-07_V7results_closed_adults_BASELINE.csv", "closedA"),
                 SA_readclean("Outputs/2023-04-27_V7results_open_adults_0.001SigmaA.csv", "0.001SigmaA"),
                 SA_readclean("Outputs/2023-05-01_V7results_open_adults_0.1SigmaA.csv", "0.1SigmaA"),
                 SA_readclean("Outputs/2023-05-01_V7results_open_adults_0.5M.csv", "0.5M"),
                 SA_readclean("Outputs/2023-05-01_V7results_open_adults_1.5M.csv", "1.5M"),
                 SA_readclean("Outputs/2023-05-01_V7results_open_adults_102sigmaR.csv", "102SigmaR"),
                 SA_readclean("Outputs/2023-05-02_V7results_open_adults_104sigmaR.csv", "104SigmaR"))
  
  
  # join data & calculate
  SAJoin <- left_join(dfSAs, dfBL,  
                        by = c("Prop_gone", "Dist_Length", "Sample_method", "Yr_Dist_Start"),
                        suffix = c("SA", "BL")) %>%
             mutate(DeltaMaxt = Samp_timeMaxSA - Samp_timeMaxBL)
  
  
  # check time to high detect in the absence of disturbance
  # (compare T2HDwoDistSA to Samp_timeMaxSA)
  # add 1 year as this is measured after the 1-year-long 'disturbance'
  T2AUCnoDist <- SAJoin %>% 
    filter(Prop_gone == 0 && Dist_Length == 1 && Yr_Dist_Start == 0) %>% 
    mutate(Samp_timeMaxSA = Samp_timeMaxSA + 1) 
  
  # stats
  T2AUCnoDistStats <- T2AUCnoDist %>% group_by(scenSA) %>% 
    summarise(mu = mean(Samp_timeMaxSA), H = max(Samp_timeMaxSA), L = min(Samp_timeMaxSA)) %>% 
    arrange(L)

  # scenario labels
  scenLabs <- c("BL" = "BASELINE",
                "BLjuvs" = "Recruits impacted by disturbance",
                "closedA" = "Closed population",
                "0.001SigmaA" = "Decreased adult variability",
                "0.1SigmaA" = "Increased adult variability",
                "0.5M" = "Decreased adult mortality",
                "1.5M" = "Increased adult mortality",
                "102SigmaR" = "Decreased recruit variability",
                "104SigmaR" = "Increased recruit variability")

  
  # scatter plot: time to HD in the absence of disturbance (this is main text figure)
  ggplot(T2AUCnoDist %>% filter(scenSA != "BLjuvs"), 
         aes(y = Samp_timeMaxSA, 
             x = factor(scenSA, levels = T2AUCnoDistStats$scenSA),
             colour = Sample_method)) +
    geom_point(size = 4, alpha = 0.8, position = position_dodge(0.1)) +
    scale_colour_manual(values = c("#1C75BC","#DA8427","#39B54A"), name = 'Sample method') +
    scale_y_continuous(name = 'Years to high detectability without a disturbance',
                       breaks = pretty_breaks()) +
    scale_x_discrete(name = "Scenario", labels = scenLabs) +
    theme_light() +
    coord_flip()
  

  
 # plot distribution of times to HD
    ggplot(SAJoin %>% 
             mutate(scenSA = factor(scenSA, levels = T2AUCnoDistStats$scenSA))) +
      geom_bar(aes(x = Samp_timeMaxSA, fill = Sample_method)) + 
      # geom_bar(aes(x = Samp_timeMaxSA, fill = Sample_method, colour = BAdistSA)) +
      facet_wrap(~scenSA, ncol = 1,
                 labeller = labeller(scenSA = scenLabs)) +
      scale_fill_manual(values = c("#1C75BC","#DA8427","#39B54A"), name = 'Sample method') +
      scale_x_continuous(name = 'Years to high detectability following disturbance', 
                         breaks = breaks_pretty(n = 15)) +
       # scale_y_discrete(name = 'Count of scenarions') +
      theme_light() +
      theme(strip.background = element_blank(),
            strip.text = element_text(colour = 'black', hjust = 0, face = "bold"))
  
  
    ggplot(SAJoin %>%
             mutate(scenSA = factor(scenSA, levels = T2AUCnoDistStats$scenSA)) %>%
             filter(scenSA != "BL")) +
      geom_bar(aes(x = DeltaMaxt, fill = Sample_method)) +
      facet_grid(rows = vars(scenSA)) +
      scale_fill_manual(values = c("#1C75BC","#DA8427","#39B54A"), name = 'Sample method') +
      scale_x_continuous(name = 'delta time to high detectability',
                         breaks = breaks_pretty(n = 15)) +
      scale_y_discrete(name = "Scenario") +
      theme_light()
  
  
  
  
  
  
  