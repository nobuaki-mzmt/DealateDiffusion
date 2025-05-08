# Data formatting for futher analysis


#------------------------------------------------------------------------------#
# This file is for preprocess all data for statistical analysis and plot
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  
  library(dplyr)

  library(lme4)
  library(car)
  library(fitdistrplus)
  library(nlme)
  
  library(CircStats)
  library(circular)
  
  library(ggplot2)
  library(patchwork)
  library(viridis)
  
  rlaplace <- function(n, location = 0, scale = 1) {
    return(location - scale * sign(runif(n) - 0.5) * log(1 - 2 * abs(runif(n) - 0.5)))
  }
  dlaplace <- function(x, location = 0, scale = 1) {
    return(1/(2 * scale) * exp(-abs(x - location) / scale))
  }
  plaplace <- function(x, location = 0, scale = 1) {
    return(0.5 * (1 + sign(x - location) * (1 - exp(-abs(x - location) / scale))))
  }
  
  # data summarize
  {
    load("data_fmt/df_expand_tandem.rda")
    df_tandem <- df
    df_summary1 <- df_tandem[df_tandem$tandem,] %>%
      group_by(x = round(step, 1)) %>%
      summarize(mean_turn = mean((turn), na.rm=T), 
                sd_turn = sd((turn), na.rm=T), 
                mean_acc = mean((acc), na.rm=T), 
                sd_acc = sd((acc), na.rm=T), 
                treat = "leader",
                .groups = "drop")
    
    load("data_fmt/df_expand_solo.rda")
    df_solo <- df
    df_summary2 <- df_solo %>%
      group_by(x = round(step, 1)) %>%
      summarize(mean_turn = mean((turn), na.rm=T), 
                sd_turn = sd((turn), na.rm=T), 
                mean_acc = mean((acc), na.rm=T), 
                sd_acc = sd((acc), na.rm=T), 
                treat = "solo",
                .groups = "drop")
    
    load("data_fmt/df_antam.rda")
    df_antam <- df
    
    df_summary3 <- df_antam %>%
      group_by(x = round(step, 1)) %>%
      summarize(mean_turn = mean((turn), na.rm=T), 
                sd_turn = sd((turn), na.rm=T), 
                mean_acc = mean((acc), na.rm=T), 
                sd_acc = sd((acc), na.rm=T), 
                treat = "antam_solo",
                .groups = "drop")
    
    df_summary <- rbind(df_summary1, df_summary2, df_summary3)
    rm(df_summary1, df_summary2, df_summary3)
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# movement parameters
#------------------------------------------------------------------------------#
{
  # speed - acceleration
  {
    # plot
    {
      acc_plot <- function(df_all, df_sum, ANTAM = FALSE, pcol = "red"){
        param <- list()
        
        if(ANTAM){ 
          r <- lmer(acc ~ step + (1|video), data = subset(df_all, step < 10))
        } else {
          r <- lmer(acc ~ step + (1|colony) + (1|colony:video), data = subset(df_all, step < 10))
        }
        intercept <- summary(r)$coefficient[1]
        slope <- summary(r)$coefficient[2]
        
        p_main <- ggplot(df_sum, aes(x, mean_acc)) +
          geom_line(color = pcol, linewidth = 1) +
          geom_ribbon(aes(ymin = mean_acc - sd_acc, 
                          ymax = mean_acc + sd_acc),
                      alpha = 0.2, fill = pcol) +
          #geom_abline(slope = slope, intercept = intercept) +
          geom_segment(data = data.frame(x = 0, xend = 10,
                                         y = intercept + slope * 0,
                                         yend = intercept + slope * 10),
                       aes(x=x, y=y, xend=xend, yend=yend),
                       color = "#2F2F2F", linetype = 3
                       ) +
          geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
          labs(y = "Change of Speed (mm / 0.2 sec)", x = "Speed (mm / 0.2 sec)") +
          theme_classic()+
          #coord_cartesian(xlim = c(0,10), ylim=c(-4,4))+
          scale_x_continuous(breaks = c(0, 5, 10), limits = c(0, 10))+
          scale_y_continuous(breaks = c(-3, 0, 3), limits = c(-4, 4))
        
        # histogram of step length
        p_step <- ggplot(df_all, aes(x = step)) +
          geom_histogram(binwidth = 0.2, fill = "#333333", color="black")+
          coord_cartesian(xlim = c(0, 10)) + 
          ylab("Count")+
          theme_classic()+
          theme(axis.text = element_text(size = 7),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                aspect.ratio = 1/5,
                panel.background = element_blank())
        hist_data <- ggplot_build(p_step)$data[[1]]
        ymax <- max(hist_data$count)
        ybreaks <- c(0, ceiling(ymax / 10000) * 10000)
        p_step <- p_step + 
          scale_y_continuous(breaks = ybreaks, limits = ybreaks)
        
        fit_acc_laplace <- fitdist(as.numeric(na.omit(df_all$acc)), 
                                   "laplace", start = list(location = 0, scale = 1))
        
        fit_laplace <- data.frame(x = seq(-5,5,0.01),
                                  y = dlaplace(seq(-5,5,0.01), 
                                               fit_acc_laplace$estimate[1],
                                               fit_acc_laplace$estimate[2]))
        
        p_acc <- ggplot(df_all, aes(y = acc)) +
          geom_histogram(binwidth = 0.2, fill = "#333333", color="black") +
          geom_path(data= fit_laplace, aes(x=y, y=x))+
          coord_cartesian(ylim = c(-4, 4)) +
          xlab("Count") + ylab("") +
          theme_classic()+
          theme(axis.text = element_text(size = 7),
                axis.title = element_blank(),
                axis.text.y = element_blank(),
                aspect.ratio = 5,
                panel.background = element_blank())
        hist_data <- ggplot_build(p_acc)$data[[1]]
        xmax <- max(hist_data$count)
        xbreaks <- c(0, ceiling(xmax / 10000) * 10000)
        p_acc <- p_acc + 
          scale_x_continuous(breaks = xbreaks, limits = xbreaks)
        
        
        print(p_step / p_main | p_acc)
        
        param$lmer <- summary(r)$coefficient
        param$fit_laplace <- fit_acc_laplace
        return(param)
      }
      
      lead_acc_param <- acc_plot(df_tandem[df_tandem$tandem,], subset(df_summary, treat == "leader"),
                                 pcol = viridis(2, option = "E")[2])
      ggsave(filename = "output/leader_acc.pdf", 
             device = cairo_pdf, family = "Arial", width = 5, height = 4)
      
      solo_acc_param <- acc_plot(df_solo, subset(df_summary, treat == "solo"),
                                 pcol = viridis(2, option = "E")[1])
      ggsave(filename = "output/solo_acc.pdf", 
             device = cairo_pdf, family = "Arial", width = 5, height = 4)
      
      if(F){
        df_antam$colony <- "Cf150704"
        df_antam$video <- df_antam$ind
        acc_plot(df_antam, subset(df_summary, treat == "antam_solo"), ANTAM = TRUE)
      }
    }
    
    # stat
    {
      # tandem vs solo
      df_combined <- bind_rows(
        df_tandem %>% filter(tandem) %>% dplyr::select(-tandem_event),
        df_solo %>% dplyr::select(-sex) %>% mutate(tandem = FALSE)
      )
      
      r <- lmer(acc ~ step * tandem + (1|colony/video), data=df_combined)
      Anova(r)
      
      r <- lmer(step ~ tandem + (1|colony/video), data=df_combined)
      Anova(r)
      
      # Model with random effect (e.g., by colony), and allowing different variance by 'tandem'
      mod_hom <- lme(acc ~ 1, 
                     random = ~ 1 | colony/video, 
                     data = df_combined,
                     na.action = na.omit)
      
      mod_het <- lme(acc ~ 1, 
                     random = ~ 1 | colony/video, 
                     weights = varIdent(form = ~ 1 | tandem),
                     data = df_combined,
                     na.action = na.omit)
      
      anova(mod_hom, mod_het)
      
      # Model with random effect (e.g., by colony), and allowing different variance by 'tandem'
      mod_hom <- lme(step ~ 1, 
                     random = ~1 | colony/video, 
                     data = df_combined,
                     na.action = na.omit)
      
      mod_het <- lme(step ~ 1, 
                     random = ~1 | colony/video, 
                     weights = varIdent(form = ~1 | tandem),
                     data = df_combined,
                     na.action = na.omit)
    
      anova(mod_hom, mod_het)
    }
  }
  
  # turning
  {
    turn_range <- c(-1, 1)
    xrange <- c(0, 10)
    # solo
    {
      turn_vals <- as.numeric(na.omit(df_solo$turn))
      fit_turn_laplace <- fitdist(turn_vals, "laplace", start = list(location = 0, scale = 1))
      fit_turn_laplace$loglik
      
      mle_wrpcauchy <- mle.wrappedcauchy(turn_vals)
      loglik_wrappedcauchy <- sum(log(dwrappedcauchy(turn_vals, as.numeric(mle_wrpcauchy$mu), mle_wrpcauchy$rho)))
      loglik_wrappedcauchy
      
      df_fit <- data.frame(x= seq(-1,1,0.01),
                 y = dlaplace(seq(-1,1,0.01), fit_turn_laplace$estimate[1], fit_turn_laplace$estimate[2]))
      
      
      p1 <- ggplot(subset(df_summary, treat == "solo"), aes(x, mean_turn)) +
        geom_line(color = viridis(2, option = "E")[1], linewidth = 1) +
        geom_ribbon(aes(ymin = mean_turn - sd_turn, 
                        ymax = mean_turn + sd_turn),
                    alpha = 0.2, fill = viridis(2, option = "E")[1]) +
        labs(y = "Turning (rad / 0.2 sec)", x = "Speed (mm / 0.2 sec)") +
        theme_classic()+
        scale_x_continuous(limits = xrange, breaks = c(0,5,10)) +
        scale_y_continuous(limits = turn_range, breaks = turn_range)
        #coord_cartesian(, ylim = turn_range)
      
      p2 <- ggplot(df_solo, aes(y = (turn), x = ..density..)) +
        geom_histogram(binwidth = 0.01, fill = "#333333", color="black") +
        geom_path(data = df_fit, aes(x = y, y = x), col = 2) +
        coord_cartesian(ylim = turn_range) +
        xlab("Count") + ylab("") +
        theme_classic()+
        theme(axis.text = element_text(size = 7),
              axis.title = element_blank(),
              axis.text.y = element_blank(),
              aspect.ratio = 5,
              panel.background = element_blank())+
        scale_y_continuous(limits = turn_range, breaks = turn_range)
      
      print(p1 | p2)
      ggsave("output/turn_solo.pdf", 
             device = cairo_pdf, family = "Arial", width = 4, height = 3)
      
      param_turn_solo <- fit_turn_laplace$estimate
    }
    
    # tandem
    {
      turn_vals <- as.numeric(na.omit(df_tandem[df_tandem$tandem,]$turn))
      fit_turn_laplace <- fitdist(turn_vals, "laplace", start = list(location = 0, scale = 1))
      fit_turn_laplace$loglik
      
      mle_wrpcauchy <- mle.wrappedcauchy(turn_vals)
      loglik_wrappedcauchy <- sum(log(dwrappedcauchy(turn_vals, as.numeric(mle_wrpcauchy$mu), mle_wrpcauchy$rho)))
      loglik_wrappedcauchy
      
      df_fit <- data.frame(x= seq(-1,1,0.01),
                           y = dlaplace(seq(-1,1,0.01), fit_turn_laplace$estimate[1], fit_turn_laplace$estimate[2]))
      
      
      p1 <- ggplot(subset(df_summary, treat == "leader"), aes(x, mean_turn)) +
        geom_line(color =  viridis(2, option = "E")[2], linewidth = 1) +
        geom_ribbon(aes(ymin = mean_turn - sd_turn, 
                        ymax = mean_turn + sd_turn),
                    alpha = 0.2, fill =  viridis(2, option = "E")[2]) +
        labs(y = "Turning (rad / 0.2 sec)", x = "Speed (mm / 0.2 sec)") +
        theme_classic()+
        scale_x_continuous(limits = xrange, breaks = c(0,5,10)) +
        scale_y_continuous(limits = turn_range, breaks = turn_range)
      #coord_cartesian(xlim = c(0,6), ylim = turn_range)
      
      p2 <- ggplot(df_tandem[df_tandem$tandem,], aes(y = (turn), x = ..density..)) +
        geom_histogram(binwidth = 0.01, fill = "#333333", color="black") +
        geom_path(data = df_fit, aes(x = y, y = x), col = 2) +
        coord_cartesian(ylim = turn_range) +
        xlab("Count") + ylab("") +
        theme_classic()+
        theme(axis.text = element_text(size = 7),
              axis.title = element_blank(),
              axis.text.y = element_blank(),
              aspect.ratio = 5,
              panel.background = element_blank())+
        scale_y_continuous(limits = turn_range, breaks = turn_range)
      #coord_cartesian(ylim = turn_range)
      print(p1 | p2)
      
      ggsave("output/turn_leader.pdf", 
             device = cairo_pdf, family = "Arial", width = 4, height = 3)
      
      param_turn_leader <- fit_turn_laplace$estimate
    }
    
    
    # Model with random effect (e.g., by colony), and allowing different variance by 'tandem'
    mod_hom <- lme(abs(turn) ~ 1, 
                   random = ~ 1 | colony/video, 
                   data = df_combined,
                   na.action = na.omit)
    
    mod_het <- lme(abs(turn) ~ 1, 
                   random = ~ 1 | colony/video, 
                   weights = varIdent(form = ~ 1 | tandem),
                   data = df_combined,
                   na.action = na.omit)
    
    anova(mod_hom, mod_het)
    
  }
  
  # parameters
  solo_param <- list(
    acc_inter = solo_acc_param$lmer[1],
    acc_slope = solo_acc_param$lmer[2],
    acc_loc   = as.numeric(solo_acc_param$fit_laplace$estimate[1]),
    acc_scale = as.numeric(solo_acc_param$fit_laplace$estimate[2]),
    turn_loc  = as.numeric(param_turn_solo[1]),
    turn_scale= as.numeric(param_turn_solo[2]))
  save(solo_param, file = "data_fmt/solo_param.rda")
  
  # traveled distance
  df_traveled_dis <- rbind(
    df_solo %>%
      group_by(video) %>%
      summarize(traveled_dis = mean(step, na.rm = T) * 5 * 60 * 15 / 1000,
                sex = unique(sex), colony = unique(colony), treat = "solo"),
    
    df_tandem[df_tandem$tandem,] %>%
      group_by(video) %>%
      summarize(traveled_dis = mean(step, na.rm = T) * 5 * 60 * 15 / 1000,
                sex = "tandem", colony = unique(colony), treat = "tandem")
  )
  
  ggplot(df_traveled_dis, aes(x = sex, y = traveled_dis, color = sex)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 2, color = "black", alpha = 0.5) +
    scale_y_continuous(limits = c(0, max(df_traveled_dis$traveled_dis)), 
                       breaks = c(0, 25)) +
    scale_color_viridis(discrete = T, 
                        option = "D", end = 0.95) +
    labs(x = "", y = "Projected traveled distance in 30 min (m)") +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1, legend.position = "none")
  
  
  ggsave("output/traveled_dis_comp.pdf", 
         device = cairo_pdf, family = "Arial", width = 4, height = 4)
  
  r <- lmer(traveled_dis ~ sex + (1|colony) , df_traveled_dis)
  Anova(r)
  mean(df_traveled_dis$traveled_dis)
  sd(df_traveled_dis$traveled_dis)
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# MSD
#------------------------------------------------------------------------------#
{
  load("data_fmt/df_msd.rda")
  
  # trajectories
  {
    prange <- 35
    ptick  <- 30
    
    ggplot(dfMSD_solo, aes(x = x/1000, y = y/1000, group = video)) +
      geom_path()+
      geom_point(data = subset(dfMSD_solo, time == 899), 
                 col = viridis(3)[2], alpha = 1, size = 1) +
      geom_point(data = subset(dfMSD_solo, time == 1798.2), 
                 col = viridis(3)[1], alpha = 1, size = 1) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = c(0.8,0.2)) +
      labs(x = "X (m)", y = "Y (m)") +
      scale_x_continuous(breaks = c(-ptick, 0, ptick), limits = c(-prange, prange)) +
      scale_y_continuous(breaks = c(-ptick, 0, ptick), limits = c(-prange, prange))
    ggsave("output/solo_trajectories_diffusion.pdf", 
           device = cairo_pdf, family = "Arial", width = 5, height = 3.5)
    
    
    ggplot(subset(dfMSD_tandem, time <= 1798.2),
           aes(x = x/1000, y = y/1000, group = tandem_event)) +
      geom_path()+
      geom_point(data = subset(dfMSD_tandem, time == 899), 
                 col = viridis(3)[2], alpha = 1, size = 1) +
      geom_point(data = subset(dfMSD_tandem, time == 1798.2), 
                 col = viridis(3)[1], alpha = 1, size = 1) +
      theme_classic() +
      theme(legend.position = c(0.8,0.2), aspect.ratio = 1) +
      labs(x = "X (m)", y = "Y (m)") +
      scale_x_continuous(breaks = c(-ptick, 0, ptick), limits = c(-prange, prange)) +
      scale_y_continuous(breaks = c(-ptick, 0, ptick), limits = c(-prange, prange))
    ggsave("output/tandem_trajectories_diffusion.pdf", 
           device = cairo_pdf, family = "Arial", width = 5, height = 3.5)
    
    ggplot(dfMSD_ANTAM, aes(x = x/1000, y = y/1000, group = ind)) +
      geom_path()+
      geom_point(data = subset(dfMSD_ANTAM, time == 899), 
                 col = viridis(3)[2], alpha = 1, size = 1) +
      geom_point(data = subset(dfMSD_ANTAM, time == 1798.2), 
                 col = viridis(3)[1], alpha = 1, size = 1) +
      theme_classic() +
      theme(aspect.ratio = 2/3, legend.position = c(0.8,0.2)) +
      labs(x = "X (m)", y = "Y (m)") +
      coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(-5, 5))
    ggsave("output/antam_trajectories_diffusion.pdf", 
           device = cairo_pdf, family = "Arial", width = 5, height = 3.5)
  }
  
  # MSD
  {
    msd_fit_ref <- c(0.2,0.4,0.6,0.8,1,2,4,6,8,10,20,40,60,80,100,200,400,600,800,1000)
  
    # solo
    df_msd_fit <- subset(dfMSD_solo, time %in% msd_fit_ref)
    r <- lmer(log10(msd) ~ log10(time) + (log10(time)|video), data = df_msd_fit)
    summary(r)
    
    ggplot(dfMSD_solo, aes(x = time, y = msd))+
      geom_path(aes(group = video), alpha = .4) +
      scale_color_viridis(discrete = T, end = .5, option = "D") +
      scale_x_log10(breaks = c(1, 10, 100, 1000)) +
      scale_y_log10(breaks = c(100, 10000, 1000000, 100000000), labels = c("10^2", "10^4", "10^6", "10^8")) +
      #coord_fixed() +
      geom_abline(slope = 1, intercept = summary(r)$coefficient[1], col = 1, linetype = 2)+
      geom_abline(slope = 2, intercept = summary(r)$coefficient[1], col = 1, linetype = 2) +
      geom_abline(slope = summary(r)$coefficient[2], intercept = summary(r)$coefficient[1], 
                  col = "firebrick", linewidth = 1) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = c(0.8,0.2)) +
      labs(x = "Tau (sec)", y = "MSD (m2)")
    ggsave("output/msd_solo.pdf", 
           device = cairo_pdf, family = "Arial", width = 3.5, height = 3.5)
    
    # tandem
    df_msd_fit <- subset(dfMSD_tandem, time %in% msd_fit_ref)
    r <- lmer(log10(msd) ~ log10(time) + (log10(time)|tandem_event), data = df_msd_fit)
    summary(r)
    
    ggplot(dfMSD_tandem, aes(x = time, y = msd))+
      geom_path(aes(group = tandem_event), alpha = .4) +
      scale_color_viridis(discrete = T, end = .5, option = "D") +
      scale_x_log10(breaks = c(1, 10, 100, 1000)) +
      scale_y_log10(breaks = c(100, 10000, 1000000, 100000000), labels = c("10^2", "10^4", "10^6", "10^8")) +
      #coord_fixed() +
      geom_abline(slope = 1, intercept = summary(r)$coefficient[1], col = 1, linetype = 2)+
      geom_abline(slope = 2, intercept = summary(r)$coefficient[1], col = 1, linetype = 2) +
      geom_abline(slope = summary(r)$coefficient[2], intercept = summary(r)$coefficient[1], 
                  col = "firebrick", linewidth = 1) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = c(0.8,0.2)) +
      labs(x = "Tau (sec)", y = "MSD (m2)")
    ggsave("output/msd_tandem.pdf", 
           device = cairo_pdf, family = "Arial", width = 3.5, height = 3.5)
    
    
    # ANMTA
    df_msd_fit <- subset(dfMSD_ANTAM, time %in% msd_fit_ref)
    r <- lmer(log10(msd) ~ log10(time) + (log10(time)|ind), data = df_msd_fit)
    
    ggplot(dfMSD_ANTAM, aes(x = time, y = msd))+
      geom_path(aes(group = ind), alpha = .4) +
      scale_color_viridis(discrete = T, end = .5, option = "D") +
      scale_x_log10(breaks = c(1, 10, 100, 1000)) +
      scale_y_log10(breaks = c(100, 10000, 1000000, 100000000), labels = c("10^2", "10^4", "10^6", "10^8")) +
      #coord_fixed() +
      geom_abline(slope = 1, intercept = summary(r)$coefficient[1], col = 1, linetype = 2)+
      geom_abline(slope = 2, intercept = summary(r)$coefficient[1], col = 1, linetype = 2) +
      geom_abline(slope = summary(r)$coefficient[2], intercept = summary(r)$coefficient[1], 
                  col = "firebrick", linewidth = 1) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = c(0.8,0.2)) +
      labs(x = "Tau (sec)", y = "MSD (m2)")
    ggsave("output/msd_antam.pdf", 
           device = cairo_pdf, family = "Arial", width = 3.5, height = 3.5)
  }
  
  # distribution
  {
    df_termite_dist <- NULL
    for(i_v in unique(dfMSD_solo$video)){
      df_temp <- subset(dfMSD_solo, video == i_v)
      df_temp$x <- df_temp$x - df_temp$x[1]
      df_temp$y <- df_temp$y - df_temp$y[1]
      df_temp <- subset(df_temp, time %in% c(899, 1798.2))
      
      df_termite_dist <- rbind(df_termite_dist, 
        data.frame(
          time_point = c(15,30),
          distance = sqrt(df_temp$x^2 + df_temp$y^2),
          treat = "solo",
          colony = df_temp$colony[1],
          video = df_temp$video[1]
        )
      )
    }
    for(i_v in unique(dfMSD_tandem$tandem_event)){
      df_temp <- subset(dfMSD_tandem, tandem_event == i_v)
      df_temp$time <- df_temp$time - df_temp$time[1]
      df_temp$x <- df_temp$x - df_temp$x[1]
      df_temp$y <- df_temp$y - df_temp$y[1]
      df_temp <- subset(df_temp, time %in% c(899, 1798.2))
      time_point = c(15,30)
      if(dim(df_temp)[1] == 0){next}
      if(dim(df_temp)[1] == 1){time_point = 15}
    
      df_termite_dist <- rbind(df_termite_dist, 
                               data.frame(
                                 time_point,
                                 distance = sqrt(df_temp$x^2 + df_temp$y^2),
                                 treat = "tandem",
                                 colony = df_temp$colony[1],
                                 video = df_temp$video[1]
                               )
      )
    }
    for(i_v in unique(dfMSD_ANTAM$ind)){
      df_temp <- subset(dfMSD_ANTAM, ind == i_v)
      df_temp$time <- df_temp$time - df_temp$time[1]
      df_temp$x <- df_temp$x - df_temp$x[1]
      df_temp$y <- df_temp$y - df_temp$y[1]
      df_temp <- subset(df_temp, time %in% c(899, 1798.2))
      time_point = c(15,30)
      if(dim(df_temp)[1] == 0){next}  
      if(dim(df_temp)[1] == 1){time_point = 15}
      df_termite_dist <- rbind(df_termite_dist, 
                               data.frame(
                                 time_point,
                                 distance = sqrt(df_temp$x^2 + df_temp$y^2),
                                 treat = "antam",
                                 colony = "antam",
                                 video = df_temp$ind[1]
                               )
      )
    }
    
    ggplot(subset(df_termite_dist, time_point == 15 & treat != "antam"), 
           aes(x = distance/1000, fill = as.factor(treat), group = treat))+
      geom_density(alpha = 0.25, adjust = 2) +
      geom_point(aes(y = ifelse(treat == "solo", -0.001, 0.001)),
                 size = 2, shape = 21) +
      #geom_histogram(alpha = 0.25, position = "identity") +
      scale_fill_viridis(discrete = T, option = "E") +
      scale_color_viridis(discrete = T, option = "E") +
      theme_classic() +
      theme(legend.position = "none", aspect.ratio = 1) +
      labs(x = "Dispersed distance (m)", y = "Density") +
      scale_y_continuous(breaks = c(0, 0.1), labels = c(0,0.1), limits = c(-0.001,0.1)) +
      scale_x_continuous(breaks = c(0, 10, 20), limits = c(0,20))
    
    ggsave("output/dispersed_dis.pdf", 
           device = cairo_pdf, family = "Arial", width = 3, height = 3)
    
    r <- lmer(distance ~ treat + (1|colony) + (1|colony:video), data= subset(df_termite_dist, time_point == 15 &treat != "antam"))
    Anova(r)
    subset(df_termite_dist, time_point == 15 &treat != "antam") %>%
      group_by(treat) %>%
      summarise(
        mean = mean(distance),
        sd = sd(distance)
      )
    
    max( subset(df_termite_dist, time_point == 15 &treat != "antam")$distance)/1000
  }
}
#------------------------------------------------------------------------------#