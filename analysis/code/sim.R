## Data analysis on tired termite project

#------------------------------------------------------------------------------#
# simulations.R
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  require(Rcpp)
  library(dplyr)
  
  library(ggplot2)
  library(viridis)
  
  load("data_fmt/solo_param.rda")
  
  Rcpp::sourceCpp("code/sim.cpp")
  
  rlaplace <- function(n, location = 0, scale = 1) {
    return(location - scale * sign(runif(n) - 0.5) * log(1 - 2 * abs(runif(n) - 0.5)))
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Encounter simulation
#------------------------------------------------------------------------------#
{
  # parameters
  {
    all_L     <- 30000
    detection <- 10
    end_time  <- 5 * 1800
    Light_L_list <- c(3000, all_L)
    num_ind_list <- c(100, 1000, 10000)
    iter <- 100
  }
  
  # simulations
  res_list <- list()
  res_pop_list <- list()
  i <- 1
  for(Light_L in Light_L_list){
    for(num_ind in num_ind_list){
      for(i_iter in 1:iter){
        print(Sys.time())
        print(paste(Light_L, num_ind, i_iter))
        res <- one_simulation(all_L, Light_L, detection, end_time, num_ind, solo_param)
        res_list[[i]] <- data.frame(
          Light_L = Light_L,
          num_ind = num_ind,
          encounter = cumsum(res$encounter),
          iter = i_iter,
          time = seq(1, end_time, 1),
          id = paste(Light_L, num_ind, i_iter, sep = "_")
        )
        res_pop_list[[i]] <- res$population
        i <- i + 1
      }
    }
  }
  df_res <- bind_rows(res_list)
  save(df_res, file = "data_fmt/sim.rda")
  
  # visualization
  load("data_fmt/sim.rda")
  
  df_res$prop <- df_res$encounter/df_res$num_ind*2
  
  df_summary <- df_res %>%
    group_by(num_ind, Light_L, time) %>%
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      sd_prop = sd(prop, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_summary$Light_L <- factor(df_summary$Light_L)
  
 ggplot(df_summary, aes(x = time/5/60, y = mean_prop, color = Light_L, fill = Light_L)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop), 
                alpha = 0.2, linetype = 0) +
    scale_color_viridis(discrete = TRUE, option = "D", end = 0.95, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "D", end = 0.95, direction = -1) +
    facet_wrap(~ num_ind, nrow = 1, ncol = 3, strip.position = "top") +
    labs(x = "Time (min)", y = "Encounter proportion") +
    scale_x_continuous(breaks = c(0,15,30)) +
    scale_y_continuous(breaks = c(0,.5,1), labels = c(0, 0.5, 1)) +
    theme_classic(base_size = 14) +
    theme(
      aspect.ratio = .85,
      legend.position = "none",
      axis.text = element_text(color = "black"),
      strip.placement = "outside",
      strip.background = element_blank()
    )
  
  ggsave("output/simulation.pdf", 
         device = cairo_pdf, family = "Arial", width = 7, height = 4)
  
  
  # example initial condition
  set.seed(123)
  num_pop <- 100
  df_initial_dis <- rbind(
    data.frame(
      x = Light_L_list[2]/2 + (runif(num_pop)-.5) * Light_L_list[1],
      y = Light_L_list[2]/2 + (runif(num_pop)-.5) * Light_L_list[1],
      treat = "light"
    ),
    data.frame(
      x = runif(num_pop) * Light_L_list[2],
      y = runif(num_pop) * Light_L_list[2],
      treat = "none"
    )
  )
  
  ggplot(df_initial_dis, aes(x, y))+
    geom_point(size = 0.05) +
    facet_wrap(~ treat, nrow = 1, ncol = 2, strip.position = "top") +
    theme_classic(base_size = 14) +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text = element_text(color = "black"),
      strip.placement = "outside",
      strip.background = element_blank()
    ) +
    labs(x = "", y = "") +
    scale_x_continuous(breaks = c(0, Light_L_list[2]), labels = c("0", "Larea"))+
    scale_y_continuous(breaks = c(0, Light_L_list[2]), labels = c("0", "Larea"))
  
  ggsave("output/sim_initial_cond.pdf", 
         device = cairo_pdf, family = "Arial", width = 7, height = 4)
  
  # SA
  # parameters
  {
    all_L_list     <- c(30000, 50000, 100000)
    Light_L_list <- c(3000, 5000, 10000, all_L)
    num_ind <- 1000
  }
  
  # simulations
  res_list <- list()
  res_pop_list <- list()
  i <- 1
  for(Light_L in Light_L_list){
    for(all_L in all_L_list){
      for(i_iter in 1:100){
        print(Sys.time())
        print(paste(all_L, Light_L, num_ind, i_iter))
        res <- one_simulation(all_L, Light_L, detection, end_time, num_ind, solo_param)
        res_list[[i]] <- data.frame(
          all_L   = all_L,
          Light_L = Light_L,
          num_ind = num_ind,
          encounter = cumsum(res$encounter),
          iter = i_iter,
          time = seq(1, end_time, 1),
          id = paste(all_L, Light_L, num_ind, i_iter, sep = "_")
        )
        res_pop_list[[i]] <- res$population
        i <- i + 1
      }
    }
  }
  df_res <- bind_rows(res_list)
  save(df_res, file = "data_fmt/sim_SA.rda")
  
  # visualization
  load("data_fmt/sim_SA.rda")
  
  df_res$prop <- df_res$encounter/df_res$num_ind*2
  
  df_summary <- df_res %>%
    group_by(num_ind, Light_L, all_L, time) %>%
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      sd_prop = sd(prop, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_summary$Light_L <- factor(df_summary$Light_L/1000)
  df_summary$all_L <- factor(df_summary$all_L/1000)
  
  ggplot(df_summary, aes(x = time/5/60, y = mean_prop, color = Light_L, fill = Light_L)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop), 
                alpha = 0.2, linetype = 0) +
    scale_color_viridis(discrete = TRUE, option = "D", end = 0.95, direction = -1) +
    scale_fill_viridis(discrete = TRUE, option = "D", end = 0.95, direction = -1) +
    facet_grid(. ~ all_L) +
    labs(x = "Time (min)", y = "Encounter proportion") +
    scale_x_continuous(breaks = c(0,15,30)) +
    scale_y_continuous(breaks = c(0,.5,1), labels = c(0, 0.5, 1)) +
    theme_classic(base_size = 14) +
    theme(
      aspect.ratio = .85,
      legend.position = "none",
      axis.text = element_text(color = "black"),
      strip.placement = "outside",
      strip.background = element_blank()
    )
  
  ggsave("output/simulation_SA.pdf", 
         device = cairo_pdf, family = "Arial", width = 7, height = 4)
  
}
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# one trajectory
#------------------------------------------------------------------------------#
{
  res <- one_simulation_plot(all_L, Light_L, detection, end_time, num_ind, solo_param)
  res$X <- (res$X + all_L/2) %% all_L
  res$Y <- (res$Y + all_L/2) %% all_L
  for(i in 1:end_time){
    par(pin=c(3,3))
    searching <- res$State[,i] == 0
    plot(res$X[searching, i], res$Y[searching, i], 
         xlim = c(0, all_L), ylim = c(0, all_L), cex = .4)
    Sys.sleep(.1)
  }
  
  {
    df <- list()
    for(id in 1:10){
      print(id)
      max_step <- 18000
      angle <- runif(1, 0, 1) * 2 * pi
      loc <- matrix(0, max_step, 2)
      speed <- runif(1, 0, 8)
      
      speed_his <- speed
      next_angle <- rlaplace(18000, location = 0, scale = solo_param$turn_scale)
      
      for(i in 2:max_step){
        speed <- speed + speed * solo_param$acc_slope + solo_param$acc_inter + rlaplace(1, 0, solo_param$acc_scale)
        speed <- pmax(speed, 0)
        speed_his <- c(speed_his, speed)
        
        angle <- angle + next_angle[i]
        
        loc[i,] <- loc[i-1,] + c(cos(angle), sin(angle)) * speed
      }
      
      msd <- c(NA, computeMSD(loc[,1], loc[,2], dim(loc)[1]-1))
      df[[id]] <- data.frame(frame = 1:max_step, time = 1:max_step/5,
                            x = loc[,1], y = loc[,2], speed_his,
                            msd, id)
    }
    
    df <- bind_rows(df)
    
    msd_fit_ref <- c(0.2,0.4,0.6,0.8,1,2,4,6,8,10,20,40,60,80,100,200,400,600,800,1000)
    df_msd_fit <- subset(df, time %in% msd_fit_ref)
    r <- lmer(log10(msd) ~ log10(time) + (log10(time)|id), data = df_msd_fit)
    ggplot(df, aes(x = time, y = msd))+
      geom_path(alpha = .4) +
      #scale_color_viridis(discrete = T, end = .5, option = "D") +
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
    
    p1 <- ggplot(df, aes(x = x, y = y, color = speed_his)) +
      geom_path()+
      scale_color_gradient(low = "blue", high = "red") +
      theme_minimal() +
      labs(color = "Value")+
      coord_fixed()
    
    p2 <- ggplot(df, aes(x = log10(frame), y = log10(x^2 + y^2))) +
      geom_point() +
      geom_abline(slope = 1, intercept = 1)+
      geom_abline(slope = 2, intercept = 1)+
      theme_classic()
    
    print(p1|p2)
  }
}
#------------------------------------------------------------------------------#
