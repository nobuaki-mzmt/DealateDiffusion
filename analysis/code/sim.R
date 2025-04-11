## Data analysis on tired termite project

#------------------------------------------------------------------------------#
# simulations.R
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  require(Rcpp)
  
  load("data_fmt/solo_param.rda")
  
  rlaplace <- function(n, location = 0, scale = 1) {
    return(location - scale * sign(runif(n) - 0.5) * log(1 - 2 * abs(runif(n) - 0.5)))
  }
}
#------------------------------------------------------------------------------#

# one trajectory
{
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
  
  df <- data.frame(frame = 1:max_step, x = loc[,1], y = loc[,2], speed_his)
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

#------------------------------------------------------------------------------#
# Simulations
#------------------------------------------------------------------------------#
run.simulation <- function(param, iter = 10000){
  
  df_sim <- NULL
  for(i_day in 0:3){
    sex_list <- c("F", "M")
    for(i_sex in 1:2){
      print(paste("Sim: sex:",sex_list[i_sex], "day:", i_day))
      
      res_vec <- rep(0, iter)
      d_searcher <- subset(df_all, day==i_day & sex==sex_list[i_sex])
      d_partner  <- subset(df_all, day==0 & sex==sex_list[3-i_sex])
      
      ids <- unique(d_searcher$id)
      sample_id_s   <- sample(ids, iter, replace = T)
      sample_time_s <- sample(1:6000, iter, replace = T)
      
      ids <- unique(d_partner$id)
      sample_id_p   <- sample(ids, iter, replace = T)
      sample_time_p <- sample(1:6000, iter, replace = T)
      
      for(i_rep in 1:iter){
        #print(paste("Sim: sex:",sex_list[i_sex], "day:", i_day, "iter:", i_rep))
        
        d_sub_s <- subset(d_searcher, id == sample_id_s[i_rep])[1:1500-1+sample_time_s[i_rep],]
        result <- randomize_trajectory(d_sub_s$x, d_sub_s$y)
        df1 <- data.frame(time = 1:1500, x=result$x, y=result$y)
        
        d_sub_p <- subset(d_partner, id == sample_id_p[i_rep])[1:1500-1+sample_time_p[i_rep],]
        result <- randomize_trajectory(d_sub_p$x, d_sub_p$y)
        df2 <- data.frame(time = 1:1500, x=result$x, y=result$y)
        
        x1 <- df1$x
        y1 <- df1$y
        x2 <- df2$x + runif(1,0,1)*L
        y2 <- df2$y + runif(1,0,1)*L
        
        res_vec[i_rep] = one_simulation(x1, y1, x2, y2, L, dis_encounter)
        #print(res_vec[i_rep])
        
      }
      df_temp <- data.frame(
        day = i_day, sex = sex_list[i_sex], encounter_time = res_vec
      )
      df_sim <- rbind(df_sim, df_temp)
    }
  }
  save(df_sim, file="data/fmt/df_sim.rda")
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# output
#------------------------------------------------------------------------------#
plot.simulations <- function(){
  load("data/fmt/df_sim.rda")
  
  df_sim$cens <- 0
  #df_sim$cens[df_sim$encounter_time==1501] <- 1
  df_sim$cens[df_sim$encounter_time == 1501] <- 1
  df_sim$encounter_time[df_sim$encounter_time == 1501] <- 3000
  
  df<-survfit(Surv(encounter_time/5)~day, type = "kaplan-meier", data=df_sim)
  ggsurvplot(fit = df, data = df_sim, fun = "event",
             pval = F, pval.method = TRUE,
             risk.table = F, conf.int = FALSE,
             ncensor.plot = FALSE, size = 1, 
             xlab="Time (sec)", 
             ylab="Encounter rate",
             palette = viridis(4),
             ggtheme = theme_classic() + theme(aspect.ratio = 1) +
               theme(strip.background = element_rect(colour="#00000000", fill="#00000000")),
             facet.by = "sex") + 
    coord_cartesian(ylim=c(0,.7), xlim = c(0,300), clip = 'on', expand=FALSE) 
  ggsave("output/simulation.pdf", width=6, height=3)  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  run.simulation(iter = 10000)
  plot.simulations()
}
#------------------------------------------------------------------------------#