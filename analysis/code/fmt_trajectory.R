# Data formatting for futher analysis


#------------------------------------------------------------------------------#
# This file is for preprocess all data for statistical analysis and plot
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  
  library(arrow)
  library(dplyr)
  
  library(CircStats)
  library(circular)
  
  library(ggplot2)
  library(patchwork)
  
  # param
  dish_size <- 145
  aframe <- 2 # frames to consider if termites approach to wall
  rad_near_wall <- 60 # how far from the center of the wall to be near wall
  
  circle_data <- data.frame(
    x = dish_size/2 * cos(seq(0, 2 * pi, length.out = 101)),
    y = dish_size/2 * sin(seq(0, 2 * pi, length.out = 101))
  )
  
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# expand dish to open space
#------------------------------------------------------------------------------#

{
  data_name <- paste0("data_fmt/solo_df.feather")
  df <- arrow::read_feather(data_name)
  df <- as.data.frame(df)
  
  # polar coordinate
  df$rad   <- sqrt(df$Center_x^2 + df$Center_y^2)
  df$theta <- atan2(df$Center_y, df$Center_x)
  
  # direction
  df$head_dir <- atan2(df$Head_y - df$Tip_y, df$Head_x - df$Tip_x)
  df$move_dir <- c(NA, atan2(diff(df$Center_y), diff(df$Center_x)))
  df$step     <- c(NA, sqrt(diff(df$Center_y)^2 + diff(df$Center_x)^2))
  
  df <- df %>% dplyr::select(Center_x, Center_y, rad, theta, head_dir, move_dir, step, video, colony, sex, frame)
  
  process_plot <- F
  
  for(i_v in unique(df$video)){
    print(i_v)
    df_temp <- df[df$video == i_v,]
    df_temp$step[1]     <- NA
    df_temp$move_dir[1] <- NA
    
    p1 <- ggplot(df_temp, aes(x = Center_x, y = Center_y))+
      geom_path(alpha = 0.5, linewidth = 0.25) +
      scale_x_continuous(limits = c(-70, 70)) + 
      scale_y_continuous(limits = c(-70, 70)) + 
      xlab("x (mm)") + ylab("y (mm)") +
      theme_classic() + 
      ggtitle(i_v) +
      theme(aspect.ratio = 1)
    
    # convert circular movements to open area
    {
      theta_diff <- diff(df_temp$theta)
      theta_diff[theta_diff > pi] <- theta_diff[theta_diff > pi] - 2 * pi
      theta_diff[theta_diff < -pi] <- theta_diff[theta_diff < -pi] + 2 * pi
      theta_diff <- c(NA, theta_diff)
      crockwise <- theta_diff < 0
      crockwise <- rle(crockwise)
      
      step_count <- 1
      for(i_w in 2:length(crockwise$lengths)){
        duration <- crockwise$length[i_w]
        wall_frames <- 1:duration + step_count
        moved_dis <- sum(df_temp[wall_frames,]$step)
        
        if(moved_dis > 10){
          
          if(process_plot){
            p1 <- ggplot(data = df_temp[wall_frames,], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
              geom_path(color = 1, alpha = 0.5) +
              geom_point(data = df_temp[wall_frames[1],], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
              xlim(-70,70) + ylim(-70,70) + 
              coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(sd(df_temp[wall_frames, ]$f_r), 2)))
          }
          
          rotation <- sum(theta_diff[wall_frames]) 
          rotation_per_frame <- (df_temp[wall_frames,]$step / sum(df_temp[wall_frames,]$step)) * rotation
          
          for(i in 1:length(wall_frames[-1])){
            df_temp2 <- df_temp
            offset_pos <- as.numeric(df_temp[wall_frames[-1][i] -1, 1:2])
            r_frames <- wall_frames[-1][i]:nrow(df_temp2)
            
            df_temp2[r_frames, 1] <- df_temp2[r_frames, 1] - offset_pos[1]
            df_temp2[r_frames, 2] <- df_temp2[r_frames, 2] - offset_pos[2]
            
            df_temp[r_frames, 1] <- cos(rotation_per_frame[i])*df_temp2[r_frames, 1] + sin(rotation_per_frame[i])*df_temp2[r_frames, 2] 
            df_temp[r_frames, 2] <- cos(rotation_per_frame[i])*df_temp2[r_frames, 2] - sin(rotation_per_frame[i])*df_temp2[r_frames, 1] 
            
            df_temp[r_frames, 1] <- df_temp[r_frames, 1] + offset_pos[1]
            df_temp[r_frames, 2] <- df_temp[r_frames, 2] + offset_pos[2]
            
            df_temp$f_dir[r_frames] <- df_temp$f_dir[r_frames] - rotation_per_frame[i]
            df_temp$f_move_dir[r_frames] <- df_temp$f_move_dir[r_frames] - rotation_per_frame[i]
          }
          
          if(process_plot){
            p2 <- ggplot(data = df_temp[wall_frames,], aes(x = fCenter_x, y = fCenter_y))+
              geom_path(color = 1, alpha = 0.5) +
              geom_point(data = df_temp[wall_frames[1],], aes(x = fCenter_x, y = fCenter_y))+
              coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(rotation))) #+ xlim(0,500) + ylim(0,500)
            ggsave(plot = p1+p2, filename = paste0("output/", i_w, ".png"))
          }
          
          
        } 
        
        step_count <- step_count + duration
      }
    }
    
    # wall bounded turn
    # look for it, based on the mismatch between head_dir and move_dir and fix it
    # coding by taking "Copfor_172_170627-10" frame 1:10 as an example
    {
      df_temp$head_dir <- atan2(sin(df_temp$head_dir), cos(df_temp$head_dir))
      df_temp$move_dir <- atan2(sin(df_temp$move_dir), cos(df_temp$move_dir))
      angle_mismatch <- df_temp$head_dir - df_temp$move_dir
      angle_mismatch <- atan2(sin(angle_mismatch), cos(angle_mismatch))
      mismatch_frame <- which(abs(angle_mismatch) > 1)
      for(i in 1:length(mismatch_frame)){
        if(mismatch_frame[i] == 1){next;}
        
        df_temp2 = df_temp1 <- df_temp
        
        # rotate the mismatch frame
        offset_pos <- as.numeric(df_temp[mismatch_frame[i]-1, 1:2])
        rotation <- angle_mismatch[mismatch_frame[i]]
        
        df_temp2[mismatch_frame[i], 1] <- df_temp2[mismatch_frame[i], 1] - offset_pos[1]
        df_temp2[mismatch_frame[i], 2] <- df_temp2[mismatch_frame[i], 2] - offset_pos[2]
        df_temp[mismatch_frame[i], 1] <- cos(rotation)*df_temp2[mismatch_frame[i], 1] + sin(rotation)*df_temp2[mismatch_frame[i], 2] 
        df_temp[mismatch_frame[i], 2] <- cos(rotation)*df_temp2[mismatch_frame[i], 2] - sin(rotation)*df_temp2[mismatch_frame[i], 1] 
        df_temp[mismatch_frame[i], 1] <- df_temp[mismatch_frame[i], 1] + offset_pos[1]
        df_temp[mismatch_frame[i], 2] <- df_temp[mismatch_frame[i], 2] + offset_pos[2]
        df_temp$move_dir[mismatch_frame[i]] <- df_temp$move_dir[mismatch_frame[i]] + rotation
        
        # add rest of the frame as the rotated place as an offset
        offset_pos <- as.numeric(df_temp1[mismatch_frame[i], 1:2])
        r_frames <- (mismatch_frame[i]+1):nrow(df_temp1)
        df_temp1[r_frames, 1] <- df_temp1[r_frames, 1] - offset_pos[1]
        df_temp1[r_frames, 2] <- df_temp1[r_frames, 2] - offset_pos[2]
        df_temp[r_frames, 1] <- df_temp1[r_frames, 1] + df_temp[mismatch_frame[i], 1]
        df_temp[r_frames, 2] <- df_temp1[r_frames, 2] + df_temp[mismatch_frame[i], 2]
        
        
      }
      df_temp$move_dir <- atan2(sin(df_temp$move_dir), cos(df_temp$move_dir))
      
      # then detect big turns near the edge and fix it
      turn_angle <- c(0, diff(df_temp$move_dir))
      turn_angle <- atan2(sin(turn_angle), cos(turn_angle))
      
      turn_near_wall <- which(df_temp$rad > rad_near_wall & abs(turn_angle) > 1)
      turn_near_wall <- turn_near_wall[turn_near_wall > aframe & turn_near_wall < dim(df_temp)[1]-aframe]
      wall_bounded_turn <- NULL
      for(i in 1:length(turn_near_wall)){
        if( mean(df_temp[(-aframe:-1)+turn_near_wall[i], ]$rad) < mean(df_temp[(1:aframe)+turn_near_wall[i], ]$rad) ){
          wall_bounded_turn <- c(wall_bounded_turn, turn_near_wall[i])
        }
      }
      
      for(i in 1:length(wall_bounded_turn)){
        df_temp2 <- df_temp
        offset_pos <- as.numeric(df_temp[wall_bounded_turn[i]-1, 1:2])
        r_frames <- wall_bounded_turn[i]:nrow(df_temp2)
        rotation <- turn_angle[wall_bounded_turn[i]]
        
        df_temp2[r_frames, 1] <- df_temp2[r_frames, 1] - offset_pos[1]
        df_temp2[r_frames, 2] <- df_temp2[r_frames, 2] - offset_pos[2]
        
        df_temp[r_frames, 1] <- cos(rotation)*df_temp2[r_frames, 1] + sin(rotation)*df_temp2[r_frames, 2] 
        df_temp[r_frames, 2] <- cos(rotation)*df_temp2[r_frames, 2] - sin(rotation)*df_temp2[r_frames, 1] 
        
        df_temp[r_frames, 1] <- df_temp[r_frames, 1] + offset_pos[1]
        df_temp[r_frames, 2] <- df_temp[r_frames, 2] + offset_pos[2]
        
        df_temp$head_dir[r_frames] <- df_temp$head_dir[r_frames] - rotation
        df_temp$move_dir[r_frames] <- df_temp$move_dir[r_frames] - rotation
        
      }
    }
    
    # plot expanded trajectories
    x_range <- range(df_temp$Center_x)
    y_range <- range(df_temp$Center_y)
    range_size <- max(diff(x_range), diff(y_range))
    center_xy <- c(mean(x_range), mean(y_range))
    x_limits <- c(center_xy[1] - range_size / 2 - dish_size, center_xy[1] + range_size / 2 + dish_size)
    y_limits <- c(center_xy[2] - range_size / 2 - dish_size, center_xy[2] + range_size / 2 + dish_size)
    
    p2 <- ggplot(df_temp, aes(x = Center_x, y = Center_y)) +
      geom_path(alpha = 0.5, linewidth = 0.25) +
      geom_polygon(data = circle_data, aes(x = x, y = y), fill = "red", alpha = 0.2) + 
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      coord_fixed() +
      xlab("x (mm)") +
      ylab(NULL) +
      theme_classic() +
      theme(aspect.ratio = 1)
    
    combined_plot <- p1 + p2 + plot_layout(guides = "collect") & theme(plot.margin = margin(.5, .5, .5, .5))
    ggsave(plot = combined_plot, width = 5.5, height = 3,
           filename = paste0("output/trajectories/", i_v, ".png"))
    
    df[df$video == i_v,] <- df_temp
  }
  
  save(df, file = "data_fmt/df_expand_solo.rda")
}



data_name <- paste0("data_fmt/tandem_df.feather")
df <- arrow::read_feather(data_name)
df <- as.data.frame(df)

df[,3:14] <- df[,3:14] - 70

# polar coordinate
df$f_r <- sqrt(df$fCenter_x^2 + df$fCenter_y^2)
df$theta <- atan2(df$fCenter_y, df$fCenter_x)

df$f_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
df$f_move_dir <- c(NA, atan2(diff(df$fCenter_y), diff(df$fCenter_x)))

df <- df %>% dplyr::select(fCenter_x, fCenter_y, f_r, theta, f_dir, f_move_dir, video, colony, frame)
df$f_step <- c(NA, sqrt(diff(df$fCenter_y)^2 + diff(df$fCenter_x)^2))

process_plot <- F

for(i_v in unique(df$video)){
  print(i_v)
  #i_v = "Copfor_172_170627-05"
  #i_v = "Copfor_172_170627-10"
  df_temp <- df[df$video == i_v,]
  df_temp$f_step[1] <- NA
  
  p1 <- ggplot(df_temp, aes(x = fCenter_x, y = fCenter_y))+
    #geom_point(alpha = 0.1) + 
    geom_path(alpha = 0.5) +
    xlim(-70,70)+ylim(-70,70)+theme(aspect.ratio = 1)
  
  theta_diff <- diff(df_temp$theta)
  theta_diff[theta_diff > pi] <- theta_diff[theta_diff > pi] - 2 * pi
  theta_diff[theta_diff < -pi] <- theta_diff[theta_diff < -pi] + 2 * pi
  theta_diff <- c(NA, theta_diff)
  crockwise <- theta_diff < 0
  crockwise <- rle(crockwise)
  
  step_count <- 1
  for(i_w in 2:length(crockwise$lengths)){
    duration <- crockwise$length[i_w]
    wall_frames <- 1:duration + step_count
    moved_dis <- sum(df_temp[wall_frames,]$f_step)
    
    if(moved_dis > 10){
      
      if(process_plot){
        p1 <- ggplot(data = df_temp[wall_frames,], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
          geom_path(color = 1, alpha = 0.5) +
          geom_point(data = df_temp[wall_frames[1],], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
          xlim(-70,70) + ylim(-70,70) + 
          coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(sd(df_temp[wall_frames, ]$f_r), 2))) #+ xlim(0,500) + ylim(0,500)
      }
      
      rotation <- sum(theta_diff[wall_frames]) 
      rotation_per_frame <- (df_temp[wall_frames,]$f_step / sum(df_temp[wall_frames,]$f_step)) * rotation
      
      for(i in 1:length(wall_frames[-1])){
        df_temp2 <- df_temp
        offset_pos <- as.numeric(df_temp[wall_frames[-1][i] -1, 1:2])
        r_frames <- wall_frames[-1][i]:nrow(df_temp2)
        
        df_temp2[r_frames, 1] <- df_temp2[r_frames, 1] - offset_pos[1]
        df_temp2[r_frames, 2] <- df_temp2[r_frames, 2] - offset_pos[2]
        
        df_temp[r_frames, 1] <- cos(rotation_per_frame[i])*df_temp2[r_frames, 1] + sin(rotation_per_frame[i])*df_temp2[r_frames, 2] 
        df_temp[r_frames, 2] <- cos(rotation_per_frame[i])*df_temp2[r_frames, 2] - sin(rotation_per_frame[i])*df_temp2[r_frames, 1] 
        
        df_temp[r_frames, 1] <- df_temp[r_frames, 1] + offset_pos[1]
        df_temp[r_frames, 2] <- df_temp[r_frames, 2] + offset_pos[2]
        
        df_temp$f_dir[r_frames] <- df_temp$f_dir[r_frames] - rotation_per_frame[i]
        df_temp$f_move_dir[r_frames] <- df_temp$f_move_dir[r_frames] - rotation_per_frame[i]
      }
      
      if(process_plot){
        p2 <- ggplot(data = df_temp[wall_frames,], aes(x = fCenter_x, y = fCenter_y))+
          geom_path(color = 1, alpha = 0.5) +
          geom_point(data = df_temp[wall_frames[1],], aes(x = fCenter_x, y = fCenter_y))+
          coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(rotation))) #+ xlim(0,500) + ylim(0,500)
        ggsave(plot = p1+p2, filename = paste0("output/", i_w, ".png"))
      }
      
      
    } 
    
    step_count <- step_count + duration
  }
  
  
  # wall bounded turn
  # look for it, based on the mismatch between f_dir and f_move_dir and fix it
  # coding by taking "Copfor_172_170627-10" frame 1:10 as an example
  df_temp$f_dir <- atan2(sin(df_temp$f_dir), cos(df_temp$f_dir))
  df_temp$f_move_dir <- atan2(sin(df_temp$f_move_dir), cos(df_temp$f_move_dir))
  angle_mismatch <- df_temp$f_dir - df_temp$f_move_dir
  angle_mismatch <- atan2(sin(angle_mismatch), cos(angle_mismatch))
  mismatch_frame <- which(abs(angle_mismatch) > 1)
  for(i in 1:length(mismatch_frame)){
    if(i == 1){next;}
    
    df_temp2 = df_temp1 <- df_temp
    
    # rotate the mismatch frame
    offset_pos <- as.numeric(df_temp[mismatch_frame[i]-1, 1:2])
    rotation <- angle_mismatch[mismatch_frame[i]]
    
    df_temp2[mismatch_frame[i], 1] <- df_temp2[mismatch_frame[i], 1] - offset_pos[1]
    df_temp2[mismatch_frame[i], 2] <- df_temp2[mismatch_frame[i], 2] - offset_pos[2]
    df_temp[mismatch_frame[i], 1] <- cos(rotation)*df_temp2[mismatch_frame[i], 1] + sin(rotation)*df_temp2[mismatch_frame[i], 2] 
    df_temp[mismatch_frame[i], 2] <- cos(rotation)*df_temp2[mismatch_frame[i], 2] - sin(rotation)*df_temp2[mismatch_frame[i], 1] 
    df_temp[mismatch_frame[i], 1] <- df_temp[mismatch_frame[i], 1] + offset_pos[1]
    df_temp[mismatch_frame[i], 2] <- df_temp[mismatch_frame[i], 2] + offset_pos[2]
    df_temp$f_move_dir[mismatch_frame[i]] <- df_temp$f_move_dir[mismatch_frame[i]] + rotation
    
    # add rest of the frame as the rotated place as an offset
    offset_pos <- as.numeric(df_temp1[mismatch_frame[i], 1:2])
    r_frames <- (mismatch_frame[i]+1):nrow(df_temp1)
    df_temp1[r_frames, 1] <- df_temp1[r_frames, 1] - offset_pos[1]
    df_temp1[r_frames, 2] <- df_temp1[r_frames, 2] - offset_pos[2]
    df_temp[r_frames, 1] <- df_temp1[r_frames, 1] + df_temp[mismatch_frame[i], 1]
    df_temp[r_frames, 2] <- df_temp1[r_frames, 2] + df_temp[mismatch_frame[i], 2]
    
    
  }
  df_temp$f_move_dir <- atan2(sin(df_temp$f_move_dir), cos(df_temp$f_move_dir))
  
  # then detect big turns near the edge and fix it
  turn_angle <- c(0, diff(df_temp$f_move_dir))
  turn_angle <- atan2(sin(turn_angle), cos(turn_angle))
  
  turn_near_wall <- which(df_temp$f_r > 60 & abs(turn_angle) > 1)
  turn_near_wall <- turn_near_wall[turn_near_wall>5]
  wall_bounded_turn <- NULL
  for(i in 1:length(turn_near_wall)){
    if( mean(df_temp[(-5:-1)+turn_near_wall[i], ]$f_r) < mean(df_temp[(1:5)+turn_near_wall[i], ]$f_r) ){
      wall_bounded_turn <- c(wall_bounded_turn, turn_near_wall[i])
    }
  }
  
  for(i in 1:length(wall_bounded_turn)){
    df_temp2 <- df_temp
    offset_pos <- as.numeric(df_temp[wall_bounded_turn[i]-1, 1:2])
    r_frames <- wall_bounded_turn[i]:nrow(df_temp2)
    rotation <- turn_angle[wall_bounded_turn[i]]
    
    df_temp2[r_frames, 1] <- df_temp2[r_frames, 1] - offset_pos[1]
    df_temp2[r_frames, 2] <- df_temp2[r_frames, 2] - offset_pos[2]
    
    df_temp[r_frames, 1] <- cos(rotation)*df_temp2[r_frames, 1] + sin(rotation)*df_temp2[r_frames, 2] 
    df_temp[r_frames, 2] <- cos(rotation)*df_temp2[r_frames, 2] - sin(rotation)*df_temp2[r_frames, 1] 
    
    df_temp[r_frames, 1] <- df_temp[r_frames, 1] + offset_pos[1]
    df_temp[r_frames, 2] <- df_temp[r_frames, 2] + offset_pos[2]
    
    df_temp$f_dir[r_frames] <- df_temp$f_dir[r_frames] - rotation
    df_temp$f_move_dir[r_frames] <- df_temp$f_move_dir[r_frames] - rotation
    
  }
  
  p2 <- ggplot(df_temp, aes(x = fCenter_x, y = fCenter_y))+
    #geom_point(alpha = 0.1) + 
    geom_path(alpha = 0.5) +
    coord_fixed() + ggtitle(i_v)
  
  ggsave(plot = p1+p2, filename = paste0(i_v, ".png"))
  
  df[df$video == i_v,] <- df_temp
}

save(df, file = "data_fmt/df_expand.rda")














# data prep
{
  load("data_fmt/df_expand.rda")
  
  f_step <- c(NA, sqrt( diff(df$fCenter_x)^2 + diff(df$fCenter_y)^2 ))
  f_turn <- diff(df$f_dir)
  f_turn[f_turn > pi] <- f_turn[f_turn > pi] - 2*pi
  f_turn[f_turn < -pi] <- f_turn[f_turn < -pi] + 2*pi
  f_turn <- c(NA, f_turn)
  
  f_turn[df$frame == 0] <- NA
  f_step[df$frame == 0] <- NA
  
  f_acc <- c(diff(f_step), NA)
  
  df_plot <- data.frame(
    frame = df$frame, 
    video = df$video, 
    colony = df$colony,
    f_step, f_turn, f_acc
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
              .groups = "drop")
  
}

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
    labs(y = "Turning", x = "Speed") +
    theme_minimal()+
    coord_cartesian(xlim = c(0,8), ylim=c(-.75,0.75))
  
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
  
}


turn <- na.omit(df_plot$f_turn)
plot(turn, c(diff(turn), NA))
plot(turn[-length(turn)], turn[-1])

x <- 1:10
x[-1]
x[-length(x)]


acf(na.omit(df_plot$f_turn))
acf(na.omit(df_plot$f_acc))
