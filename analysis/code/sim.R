## Data analysis on tired termite project

#------------------------------------------------------------------------------#
# simulations.R
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  require(Rcpp)
  
  library(fitdistrplus)
  rlaplace <- function(n, location = 0, scale = 1) {
    return(location - scale * sign(runif(n) - 0.5) * log(1 - 2 * abs(runif(n) - 0.5)))
  }
  dlaplace <- function(x, location = 0, scale = 1) {
    return(1/(2 * scale) * exp(-abs(x - location) / scale))
  }
  plaplace <- function(x, location = 0, scale = 1) {
    return(0.5 * (1 + sign(x - location) * (1 - exp(-abs(x - location) / scale))))
  }
  fit_acceleration_laplace <- fitdist(as.numeric(na.omit(df_plot$f_acc)), "laplace", start = list(location = 0, scale = 1))
  
  fit_turn_laplace <- fitdist(as.numeric(na.omit(df_plot$f_turn)), "laplace", start = list(location = 0, scale = 1))
  plot(density(as.numeric(na.omit(df_plot$f_turn))))
  points(seq(-1,1,0.01), dlaplace(seq(-1,1,0.01), 0, 0.12835), col =2, type ="l")
}
#------------------------------------------------------------------------------#

max_step <- 18000
angle <- runif(1, 0, 1) * 2 * pi
loc <- matrix(0, max_step, 2)
speed <- runif(1, 0, 8)

speed_his <- speed
next_angle_1 <- as.numeric(rwrappedcauchy(max_step, pi, param$angle_rho_1+0.02) - pi)
while(sum(abs(next_angle_1) > 0.8)>0){
  print(sum(abs(next_angle_1) > 0.8))
  next_angle_1[abs(next_angle_1) > 0.8] <- 
    as.numeric(rwrappedcauchy(length(next_angle_1[abs(next_angle_1) > 0.8]), pi, param$angle_rho_1+0.02) - pi)
}
next_angle_2 <- as.numeric(rwrappedcauchy(max_step, pi, param$angle_rho_2+0.02) - pi)
while(sum(abs(next_angle_2) > 0.8)>0){
  next_angle_2[abs(next_angle_2) > 0.8] <- 
    as.numeric(rwrappedcauchy(length(next_angle_2[abs(next_angle_2) > 0.8]), pi, param$angle_rho_2+0.02) - pi)
}
next_angle_3 <- as.numeric(rwrappedcauchy(max_step, pi, param$angle_rho_3+0.02) - pi)
while(sum(abs(next_angle_3) > 0.8)>0){
  next_angle_3[abs(next_angle_3) > 0.8] <- 
    as.numeric(rwrappedcauchy(length(next_angle_3[abs(next_angle_3) > 0.8]), pi, param$angle_rho_3+0.02) - pi)
}

speed_his <- speed

for(i in 2:max_step){
  speed <- speed + speed * param$slope + param$intercept + rlaplace(1, 0, fit_acceleration_laplace$estimate[2])
  speed <- pmax(speed, 0)
  speed_his <- c(speed_his, speed)
  
  #next_angle <- rwrappedcauchy(1, pi, 0.925) - pi
  next_angle <- rlaplace(1, 0, 0.12835)
  #if(speed <= 1){
  #  next_angle <- next_angle_1[i]
  #} else if (speed <= 3){
  #  next_angle <- next_angle_2[i]
  #} else{
  #  next_angle <- next_angle_3[i]
  #}
  next_angle <- pmin(next_angle, 1)
  next_angle <- pmax(next_angle, -1)
  angle <- angle + next_angle
  
  loc[i,] <- loc[i-1,] + c(cos(angle), sin(angle))*speed
  #plot(loc[i,1], loc[i,2], xlim=c(-200,200), ylim=c(-200,200))
  #Sys.sleep(0.1)
}

df <- data.frame(frame = 1:max_step, x = loc[,1], y = loc[,2], speed_his)
ggplot(df, aes(x = x, y = y, color = speed_his)) +
  geom_path()+
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(color = "Value")+
  coord_fixed()

ggplot(df, aes(x = log10(frame), y = log10(x^2 + y^2))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 1)+
  geom_abline(slope = 2, intercept = 1)

sum(speed_his)/1000


df_plot

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