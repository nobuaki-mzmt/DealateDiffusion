# Correlated random walk model with sleap
# N Mizumoto
{
  rm(list = ls())
  
  library(arrow)
  library(dplyr)
  
  library(CircStats)
  library(circular)
  
  library(ggplot2)
  
}


# data prep
{
  data_name <- paste0("data_fmt/tandem_df.feather")
  df <- arrow::read_feather(data_name)
  df <- as.data.frame(df)
  
  
  f_step_x <- c(NA, diff(df$fCenter_x))
  f_step_y <- c(NA, diff(df$fCenter_y))
  
  f_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
  m_dir <- atan2(df$mHead_y - df$mTip_y, df$mHead_x - df$mTip_x)
  
  f_step <- c(NA, sqrt( diff(df$fCenter_x)^2 + diff(df$fCenter_y)^2 ))
  m_step <- c(NA, sqrt( diff(df$mCenter_x)^2 + diff(df$mCenter_y)^2 ))
  f_turn <- diff(f_dir)
  m_turn <- diff(m_dir)
  f_turn[f_turn > pi] <- f_turn[f_turn > pi] - 2*pi
  m_turn[m_turn > pi] <- m_turn[m_turn > pi] - 2*pi
  f_turn[f_turn < -pi] <- f_turn[f_turn < -pi] + 2*pi
  m_turn[m_turn < -pi] <- m_turn[m_turn < -pi] + 2*pi
  f_turn <- c(NA, f_turn)
  m_turn <- c(NA, m_turn)
  
  f_turn[df$frame == 0] <- NA
  m_turn[df$frame == 0] <- NA
  f_step[df$frame == 0] <- NA
  m_step[df$frame == 0] <- NA
  f_step_x[df$frame == 0] <- NA
  f_step_y[df$frame == 0] <- NA
  
  #f_turn <- abs(f_turn)
  #m_turn <- abs(m_turn)
  
  f_acc <- c(diff(f_step), NA)
  m_acc <- c(diff(m_step), NA)
  
  df_plot <- data.frame(
    frame = df$frame, 
    video = df$video, 
    colony = df$colony,
    f_step, m_step, f_turn, m_turn, f_acc, m_acc, wall, f_step_x, f_step_y
  )
  
  #df_plot <- df_plot[!df_plot$wall,]
  #df_plot$f_turn [!is.na(df_plot$f_turn) & abs(df_plot$f_turn) < 0.083 ] <- 0
  #df_plot$f_turn [!is.na(df_plot$f_turn) & df_plot$f_turn > 0.083 ] <- df_plot$f_turn [!is.na(df_plot$f_turn) & df_plot$f_turn > 0.083 ] - 0.083
  #df_plot$f_turn [!is.na(df_plot$f_turn) & df_plot$f_turn < -0.083] <- df_plot$f_turn [!is.na(df_plot$f_turn) & df_plot$f_turn < -0.083] + 0.083
  
  df_summary <- df_plot %>%
    group_by(x = round(f_step, 1)) %>%
    summarize(mean_f_turn = mean((f_turn), na.rm=T), 
              sd_f_turn = sd((f_turn), na.rm=T), 
              mean_f_acc = mean((f_acc), na.rm=T), 
              sd_f_acc = sd((f_acc), na.rm=T), 
              mean_m_turn = mean((m_turn), na.rm=T), 
              sd_m_turn = sd((m_turn), na.rm=T), 
              mean_m_acc = mean((m_acc), na.rm=T), 
              sd_m_acc = sd((m_acc), na.rm=T), 
              .groups = "drop")
  
}


tapply(subset(df_plot, f_step_x > 0)$f_step_x, subset(df_plot, f_step_x > 0)$video, sum) / 140
tapply(subset(df_plot, f_step_x < 0)$f_step_x, subset(df_plot, f_step_x < 0)$video, sum)
tapply(subset(df_plot, f_step_y > 0)$f_step_y, subset(df_plot, f_step_y > 0)$video, sum) / 140
tapply(subset(df_plot, f_step_y < 0)$f_step_y, subset(df_plot, f_step_y < 0)$video, sum)

sum((tapply(subset(df_plot, f_step_x > 0)$f_step_x, subset(df_plot, f_step_x > 0)$video, sum) / 280 +
  tapply(subset(df_plot, f_step_y > 0)$f_step_y, subset(df_plot, f_step_y > 0)$video, sum) / 280)/2)

3651.857 * 2 * pi

200*2*pi / (5*60*60)

res <- NULL
for(i in seq(0.7, 1, 0.01)){
  res <- c(res,sum(abs(rwrappedcauchy(5*60*60, pi, i) - pi)))
  
}
plot(seq(0.7, 1, 0.01), res)

cbind(seq(0.7, 1, 0.01), res)

diff(res)

rwrappedcauchy(100, pi, 0.91) - pi

# plot kinetics
param <- list()
{
  # leader steps
  hist_info <- hist(df_plot$f_step, 
                    breaks = seq(min(df_plot$f_step, na.rm = TRUE), 
                                 max(df_plot$f_step, na.rm = TRUE) + 0.1, 
                                 by = 0.1), 
                    plot = F)
  hist_info$breaks[c(TRUE, hist_info$count > 100)]
  
  # histogram of female step
  ggplot(df_plot, aes(x = f_step)) +
    geom_histogram(binwidth = 0.1, fill = "#333333", color="black")+
    coord_cartesian(xlim = c(0, 8.6)) + 
    ylab("Count")+
    theme_classic()+
    theme(axis.text = element_text(size = 7),
          axis.title = element_blank(),
          #axis.text.x = element_blank(),
          aspect.ratio = 1/5,
          panel.background = element_blank())
  
  
  # move / pause 
  {
    res <- rle(f_step < 0.1)
    res$lengths <-  res$lengths[!is.na(res$values)]
    res$values <- res$values[!is.na(res$values)]
    truehist(res$lengths[res$values])
    pause <- res$lengths[res$values]
    
    library(poweRlaw)
    
    pl_model <-  displ$new(pause)
    exp_model <- disexp$new(pause)
    
    2-2*log(estimate_pars(pl_model)$value)
    2-2*log(estimate_pars(exp_model)$value)
    
    estimate_pars(exp_model)$pars
    
    # exponential -> so it is just random process
    # we do not think specific pausing behavior here
    }
  
  # speed - acceleration
  {
    library(lme4)
    library(car)
    r <- lmer(f_acc ~ f_step + (1|colony) + (1|colony:video), data = df_plot)
    intercept <- summary(r)$coefficient[1]
    slope <- summary(r)$coefficient[2]
    
    ggplot(df_summary, aes(x, mean_f_acc)) +
      geom_line(color = "red", linewidth = 1) +
      geom_ribbon(aes(ymin = mean_f_acc - sd_f_acc, 
                      ymax = mean_f_acc + sd_f_acc),
                  alpha = 0.2, fill = "red") +
      geom_abline(slope = slope, intercept = intercept) +
      labs(y = "acc", x = "Speed") +
      theme_minimal()+
      coord_cartesian(xlim = c(0,8.6), ylim=c(-2,2))
    
    plot(sd_f_acc ~ x, data = subset(df_summary, x < 8.6), xlim=c(0,8))
    
    library(ggridges)
    df_sd <- subset(df_plot, f_step < 8.6)
    df_sd$f_step <- factor(round(df_sd$f_step))
    ggplot(df_sd, aes(x = f_acc, y = f_step, group = f_step)) +
      geom_density_ridges() +
      coord_cartesian(xlim = c(-4, 4)) +
      theme_classic()+
      theme(axis.text = element_text(size = 7),
            axis.title = element_blank(),
            #axis.text.x = element_blank(),
            aspect.ratio = 1,
            panel.background = element_blank())
    
    fligner.test(f_acc ~ f_step, data = df_sd)
    sd(df_plot$f_acc, na.rm=T)
    rnorm(1000, 0, 0.8274498)
    
    truehist(df_plot$f_acc)
    points(seq(-5,5,0.01), dnorm(seq(-5,5,0.01), 0 , sd(df_plot$f_acc, na.rm=T)), type = "l")    
    
    fit_acceleration_laplace <- fitdist(as.numeric(na.omit(df_plot$f_acc)), "laplace", start = list(location = 0, scale = 1))
    fit_acceleration_laplace
    points(seq(-5,5,0.01), dlaplace(seq(-5,5,0.01), fit_acceleration_laplace$estimate[1],
                                    fit_acceleration_laplace$estimate[2]), type = "l")    
    
    rlaplace(1, 0, fit_acceleration_laplace$estimate[2])
    
    param$intercept <- intercept
    param$slope <- slope
    param$acc_sd <- sd(df_plot$f_acc, na.rm=T)
  }
  
  # turning
  ggplot(df_summary, aes(x, mean_f_turn)) +
    geom_line(color = "red", linewidth = 1) +  # Mean line
    geom_ribbon(aes(ymin = mean_f_turn - sd_f_turn, 
                    ymax = mean_f_turn + sd_f_turn),
                alpha = 0.2, fill = "red") +
    geom_ribbon(aes(ymin = mean_f_turn - 0.031, 
                    ymax = mean_f_turn + 0.031),
                alpha = 0.2, fill = "blue") +
    labs(y = "Turning", x = "Speed") +
    theme_minimal()+
    coord_cartesian(xlim = c(0,8), ylim=c(-.5,0.5))
  
  plot(sd_f_turn ~ x, data = df_summary, xlim=c(0,8))
  
  df_angle <- df_plot %>%
    mutate(angle_group = case_when(
      f_step <= 1 ~ "f_step <= 1",
      f_step > 1 & f_step <= 3 ~ "1 < f_step <= 3",
      f_step > 3 ~ "f_step > 3"
    ))
  
  ggplot(df_angle[df_angle$f_turn!=0,], aes(x = f_turn, y = angle_group)) +
    geom_density_ridges()+
    theme_classic()+
    coord_cartesian(xlim = c(-pi,pi)) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_blank(),
          #axis.text.x = element_blank(),
          aspect.ratio = 1,
          panel.background = element_blank())
  
  
  get_rho <- function(angle_data){
    mle_vonmises <- mle.vonmises(angle_data)
    mle_wrpcauchy <- mle.wrappedcauchy(angle_data)
    
    loglik_vm <- sum(dvonmises(angle_data, mle_vonmises$mu, mle_vonmises$kappa))
    loglik_wc <- sum(dwrappedcauchy(angle_data, as.numeric(mle_wrpcauchy$mu), mle_wrpcauchy$rho))
    
    plot(density(angle_data))
    points(seq(-pi,pi,0.001), circular::dwrappedcauchy(seq(-pi,pi,0.001), mu = mle_wrpcauchy$mu, rho = mle_wrpcauchy$rho), type="l", col=3)
    points(seq(-pi,pi,0.001), circular::dvonmises(seq(-pi,pi,0.001), mu = mle_vonmises$mu, kappa = mle_vonmises$kappa), type="l", col = 2)
    return(mle_wrpcauchy$rho)
  }
  
  angle_data <- na.omit(subset(df_plot, f_step < 1 & f_turn!=0)$f_turn)
  param$angle_rho_1 <- get_rho(angle_data)
  angle_data <- na.omit(subset(df_plot, f_step < 3 & f_step > 1 & f_turn!=0)$f_turn)
  param$angle_rho_2 <- get_rho(angle_data)
  angle_data <- na.omit(subset(df_plot, f_step > 3 & f_step < 8 & f_turn!=0)$f_turn)
  param$angle_rho_3 <- get_rho(angle_data)
  
  param$angle_rho <- get_rho(na.omit(df_plot$f_turn))
  
  sum(na.omit(df_plot$f_turn)) - (1825.929 * 2 * pi)
  
  (1825.929 * 2 * pi)/length(na.omit(df_plot$f_turn))
  
  optimize(ff, c(0.8,0.9999))
  ff <- function(x){ 
    abs(sum( abs(rwrappedcauchy(length(na.omit(df_plot$f_turn)), pi, x) - pi )) - (1825.929 * 2 * pi))
  }
  
  plot(density(angle_data), ylim=c(0,30))
  angle_data2 <- na.omit(df_plot$f_turn)
  #angle_data2[abs(angle_data2) < 0.032] <- 0
  angle_data2[angle_data2 > 0.0543524] <- angle_data2[angle_data2 > 0.0543524] - 0.0543524
  angle_data2[angle_data2 < -0.0543524] <- angle_data2[angle_data2 < -0.0543524] + 0.0543524
  get_rho(angle_data2)
  
  points(seq(-pi,pi,0.01), circular::dwrappedcauchy(seq(-pi,pi,0.01), mu = mle_wrpcauchy$mu, rho = 0.9922753), type="l", col=3)
  
  
}


param


