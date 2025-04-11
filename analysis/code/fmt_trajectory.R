# Data formatting for futher analysis


#------------------------------------------------------------------------------#
# This file is for preprocess all data for statistical analysis and plot
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  rm(list = ls())
  
  library(arrow)
  library(dplyr)
  library(data.table)
    
  library(ggplot2)
  library(patchwork)
  library(viridis)
  
  
  # param
  dish_size <- 145
  aframe        <- 2   # frames to consider if termites approach to wall
  rad_near_wall <- 60  # how far from the center of the wall to be near wall
  tandem_dis    <- 3.1 # ftip-mhead distance threshold based on Mizumoto and Reiter 2025 
  
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
  ### Solo
  {
    df <- arrow::read_feather("data_fmt/solo_df.feather")
    df <- as.data.frame(df)
    
    df <- df[df$frame < 53947,]

    # polar coordinate
    df$rad   <- sqrt(df$Center_x^2 + df$Center_y^2)
    df$theta <- atan2(df$Center_y, df$Center_x)
    
    # direction
    df$head_dir <- atan2(df$Head_y - df$Tip_y, df$Head_x - df$Tip_x)
    df$move_dir <- c(NA, atan2(diff(df$Center_y), diff(df$Center_x)))
    df$step     <- c(NA, sqrt(diff(df$Center_y)^2 + diff(df$Center_x)^2))
    
    #df$turn <- c(NA, diff(df$head_dir))
    #df$turn <- atan2(sin(df$turn), cos(df$turn))
    #df[abs(df$turn) > 1.5 & df$frame > 0,]
    
    
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
          
          if(duration > 1 & moved_dis > 10){
            
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
              
              df_temp$f_head_dir[r_frames] <- df_temp$f_head_dir[r_frames] - rotation_per_frame[i]
              df_temp$f_move_dir[r_frames] <- df_temp$f_move_dir[r_frames] - rotation_per_frame[i]
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
        mismatch_frame <- which(abs(angle_mismatch) > 0.5)
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
        if(length(turn_near_wall) > 0){
          for(i in 1:length(turn_near_wall)){
            if( mean(df_temp[(-aframe:-1)+turn_near_wall[i], ]$rad) < mean(df_temp[(1:aframe)+turn_near_wall[i], ]$rad) ){
              wall_bounded_turn <- c(wall_bounded_turn, turn_near_wall[i])
            }
          }
          if(length(wall_bounded_turn) > 0){
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
        }
      }
      
      # plot expanded trajectories
      x_range <- range(df_temp$Center_x, na.rm = T)
      y_range <- range(df_temp$Center_y, na.rm = T)
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
    
    df$acc  <- c(diff(df$step), NA)
    df$turn <- c(NA, diff(df$head_dir))
    df$turn <- atan2(sin(df$turn), cos(df$turn))
    df$turn[df$frame == 0] <- NA
    df_msd <- df %>% dplyr::select(frame, x = Center_x, y = Center_y, video, colony, sex)
    df <- df %>% dplyr::select(frame, head_dir, move_dir, step, acc, turn, video, colony, sex)
    
    save(df, df_msd, file = "data_fmt/df_expand_solo.rda")
  }
  
  ### Tandem
  {
    df <- arrow::read_feather("data_fmt/tandem_df.feather")
    df <- as.data.frame(df)
    
    df[,3:14] <- df[,3:14] - dish_size/2
    
    # tandem running behavior
    {
      ftoM_dis  <- sqrt( (df$fTip_x-df$mHead_x)^2 + (df$fTip_y-df$mHead_y)^2)
      tandem <- ftoM_dis < tandem_dis
      tandem_duration <- rle(tandem)
      
      # remove short separation
      count <- 1
      for(i in 1:length(tandem_duration$lengths)){
        t_duration <- tandem_duration$lengths[i]
        t_or_not   <- tandem_duration$values[i]
        if(!t_or_not && t_duration < 5){
          tandem[count:(count + t_duration - 1)] <- TRUE
        }
        count <- count + t_duration
      }
      
      # give event number
      tandem_duration <- rle(tandem)
      tandem_event <- tandem
      count <- 1
      event_count <- 1
      for(i in 1:length(tandem_duration$lengths)){
        t_duration <- tandem_duration$lengths[i]
        tandem_event[count:(count + t_duration - 1)] <- event_count
        count <- count + t_duration
        event_count <- event_count + 1
      }
      df$tandem <- tandem
      df$tandem_event <- tandem_event
    }
    
    # polar coordinate
    df$f_r <- sqrt(df$fCenter_x^2 + df$fCenter_y^2)
    df$theta <- atan2(df$fCenter_y, df$fCenter_x)
    
    df$f_head_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
    df$f_move_dir <- c(NA, atan2(diff(df$fCenter_y), diff(df$fCenter_x)))
    df$f_step <- c(NA, sqrt(diff(df$fCenter_y)^2 + diff(df$fCenter_x)^2))
    
    df <- df %>% dplyr::select(fCenter_x, fCenter_y, f_r, theta, f_head_dir, f_move_dir, f_step, tandem, tandem_event, video, colony, frame)
    
    process_plot <- F
    
    for(i_v in unique(df$video)){
      
      print(i_v)
      # i_v = "Copfor_172_170627-01"
      # i_v = "Copfor_172_170627-05"
      # i_v = "Copfor_172_170627-10"
      df_temp <- df[df$video == i_v,]
      df_temp$f_step[1] <- NA
      
      p1 <- ggplot(df_temp, aes(x = fCenter_x, y = fCenter_y, col = tandem, group = tandem_event))+
        #geom_point(alpha = 0.1) + 
        geom_path(alpha = 0.5, linewidth = 0.25) +
        scale_color_viridis(discrete = T, end = 0.8) +
        scale_x_continuous(limits = c(-70, 70)) + 
        scale_y_continuous(limits = c(-70, 70)) + 
        xlab("x (mm)") + ylab("y (mm)") +
        theme_classic() + 
        ggtitle(i_v) +
        theme(aspect.ratio = 1, legend.position = "none")
      
      # expansion
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
          moved_dis <- sum(df_temp[wall_frames,]$f_step)
          
          if(moved_dis > 10 & duration > 1){
            
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
              
              df_temp$f_head_dir[r_frames] <- df_temp$f_head_dir[r_frames] - rotation_per_frame[i]
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
        # look for it, based on the mismatch between f_head_dir and f_move_dir and fix it
        # coding by taking "Copfor_172_170627-10" frame 1:10 as an example
        df_temp$f_head_dir <- atan2(sin(df_temp$f_head_dir), cos(df_temp$f_head_dir))
        df_temp$f_move_dir <- atan2(sin(df_temp$f_move_dir), cos(df_temp$f_move_dir))
        angle_mismatch <- df_temp$f_head_dir - df_temp$f_move_dir
        angle_mismatch <- atan2(sin(angle_mismatch), cos(angle_mismatch))
        mismatch_frame <- which(abs(angle_mismatch) > 0.5)
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
        
        turn_near_wall <- which(df_temp$f_r > rad_near_wall & abs(turn_angle) > 1)
        turn_near_wall <- turn_near_wall[turn_near_wall > aframe & turn_near_wall < dim(df_temp)[1]-aframe]
        wall_bounded_turn <- NULL
        for(i in 1:length(turn_near_wall)){
          if( mean(df_temp[(-aframe:-1)+turn_near_wall[i], ]$f_r) < mean(df_temp[(1:aframe)+turn_near_wall[i], ]$f_r) ){
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
          
          df_temp$f_head_dir[r_frames] <- df_temp$f_head_dir[r_frames] - rotation
          df_temp$f_move_dir[r_frames] <- df_temp$f_move_dir[r_frames] - rotation
          
        }
      }
      
      # plot expanded trajectories
      x_range <- range(df_temp$fCenter_x)
      y_range <- range(df_temp$fCenter_y)
      range_size <- max(diff(x_range), diff(y_range))
      center_xy <- c(mean(x_range), mean(y_range))
      x_limits <- c(center_xy[1] - range_size / 2 - dish_size, center_xy[1] + range_size / 2 + dish_size)
      y_limits <- c(center_xy[2] - range_size / 2 - dish_size, center_xy[2] + range_size / 2 + dish_size)
      
      
      p2 <- ggplot() +
        geom_path(data = df_temp, aes(x = fCenter_x, y = fCenter_y, 
                                      col = tandem, group = tandem_event),
                  alpha = 0.5, linewidth = 0.25) +
        geom_path(data = df_temp[!df_temp$tandem,], 
                  aes(x = fCenter_x, y = fCenter_y, group = tandem_event),
                  alpha = 1, linewidth = 0.4, col = viridis(2)[1]) +
        #geom_point(data = df_temp[!df_temp$tandem,], 
        #           aes(x = fCenter_x, y = fCenter_y),
        #          alpha = 0.1, size = 1, shape = 16, col = viridis(2)[1]) +
        geom_polygon(data = circle_data, aes(x = x, y = y), 
                     fill = "red", alpha = 0.2) + 
        scale_color_viridis(discrete = T, end = 0.8) +
        scale_x_continuous(limits = x_limits) +
        scale_y_continuous(limits = y_limits) +
        coord_fixed() +
        xlab("x (mm)") +
        ylab(NULL) +
        theme_classic() +
        theme(aspect.ratio = 1, legend.position = "none")

      combined_plot <- p1 + p2 + plot_layout(guides = "collect") & theme(plot.margin = margin(.5, .5, .5, .5))
      ggsave(plot = combined_plot, width = 5.5, height = 3,
             filename = paste0("output/trajectories/", i_v, ".png"))
      
      df[df$video == i_v,] <- df_temp
    }
    
    
    df$acc  <- c(diff(df$f_step), NA)
    df$turn <- c(NA, diff(df$f_head_dir))
    df$turn <- atan2(sin(df$turn), cos(df$turn))
    df$turn[df$frame == 0] <- NA
    
    df_msd <- df %>% dplyr::select(frame, x = fCenter_x, y = fCenter_y, 
                                   tandem, tandem_event, video, colony)
    df <- df %>% dplyr::select(frame, head_dir = f_head_dir, move_dir = f_move_dir,
                               step = f_step, acc, turn, video, colony, tandem, tandem_event)
    
    save(df, df_msd, file = "data_fmt/df_expand_tandem.rda")
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# preprocess ANTAM data
#------------------------------------------------------------------------------#
{
  f.namesplace <- list.files("ANTAM/data_raw", pattern=".csv",full.names=T)
  f.names <- list.files("ANTAM/data_raw", pattern=".csv",full.names=F)
  
  df <- NULL
  for(i_files in 1:length(f.namesplace)){
    print(paste0(i_files, "/", length(f.names), " ", f.names[i_files]))
    
    # prep data
    d <- data.frame(fread(f.namesplace[i_files], header=F))[,c(4,1,2)]
    diff(plot(d[,1]))
    
    colnames(d) = c("time", "x", "y")
    
    #print(ggplot(d, aes(x,y))+geom_path())
  
    
    sex <- str_split(f.names[i_files], "_")[[1]][4]
    ind <- str_remove(f.names[i_files], ".csv")
    
    d[,2:3] <- d[,2:3] * 0.023725 # scale
    d <- rbind(c(0,0,0), d)
    
    # smoothing
    d$x = runmed(d$x, 5)
    d$y = runmed(d$y, 5)
    
    # interpolate
    time = seq(0, 1800, 0.2)
    x = approx(d$time, d$x, xout=time, method = "linear")$y
    y = approx(d$time, d$y, xout=time, method = "linear")$y
    
    #
    step         <- sqrt( diff(x)^2 + diff(y)^2)
    traveled_dis <- cumsum(step)
    turn         <- atan2(diff(y), diff(x))
    sqrdispl     <- x^2 + y^2
    
    df_temp <- data.frame(
      time, x, y,
      step = c(NA, step), 
      acc  = c(NA, diff(step), NA),
      turn = c(NA, turn),  
      traveled_dis = c(0, traveled_dis),
      sqrdispl, sex, ind
    )
    df <- rbind(df, df_temp)
    
    # plot trajectories
    x_range <- range(df_temp$x)
    y_range <- range(df_temp$y)
    range_size <- max(diff(x_range), diff(y_range))
    center_xy <- c(mean(x_range), mean(y_range))
    x_limits <- c(center_xy[1] - range_size / 2, center_xy[1] + range_size / 2)
    y_limits <- c(center_xy[2] - range_size / 2, center_xy[2] + range_size / 2)
    
    ggplot(df_temp, aes(x = x, y = y)) + 
      geom_path(alpha = 0.5, linewidth = 0.25) +
      geom_point(data = data.frame(x = 0, y = 0), aes(x,y), col = 2, alpha = 0.5) +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      coord_fixed() +
      xlab("x (mm)") +
      ylab(NULL) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = "none") +
      ggtitle(str_remove(ind, "mouselog_"))
    ggsave(width = 3, height = 3,
           filename = paste0("output/trajectories/antam_", str_remove(ind, "mouselog_"), ".png"))
  }
  
  save(df, file = "data_fmt/df_antam.rda")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# plot individual speed development
#------------------------------------------------------------------------------#
{
  load("data_fmt/df_expand_tandem.rda")
  for(i_v in unique(df$video)){
    df_temp <- subset(df, video == i_v)
    p1 <- ggplot(df_temp, aes(x = frame)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_path(aes(y = step), alpha = 0.4) + 
      geom_path(aes(y = acc), alpha = 0.4, col = 2) +
      geom_path(aes(y = turn), alpha = 0.4, col = 3) +
      theme_classic() +
      coord_cartesian(ylim = c(-5,10)) +
      ylab("Variable (mm or rad)") +
      ggtitle(i_v)
    ggsave(paste0("output/kinetics/", i_v, ".png"), plot = p1, width = 4, height = 3)
  }
  
  load("data_fmt/df_expand_solo.rda")
  for(i_v in unique(df$video)){
    df_temp <- subset(df, video == i_v)
    p1 <- ggplot(df_temp, aes(x = frame)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_path(aes(y = step), alpha = 0.4) + 
      geom_path(aes(y = acc), alpha = 0.4, col = 2) +
      geom_path(aes(y = turn), alpha = 0.4, col = 3) +
      theme_classic() +
      coord_cartesian(ylim = c(-5,10)) +
      ylab("Variable (mm or rad)") +
      ggtitle(i_v)
    ggsave(paste0("output/kinetics/", i_v, ".png"), plot = p1, width = 4, height = 3)
  }
  
  load("data_fmt/df_antam.rda")
  for(i_v in unique(df$ind)){
    df_temp <- subset(df, ind == i_v)
    p1 <- ggplot(df_temp, aes(x = time)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_path(aes(y = step), alpha = 0.4) + 
      geom_path(aes(y = acc), alpha = 0.4, col = 2) +
      geom_path(aes(y = turn), alpha = 0.4, col = 3) +
      theme_classic() +
      coord_cartesian(ylim = c(-5,10)) +
      ylab("Variable (mm or rad)") +
      ggtitle(i_v)
    ggsave(paste0("output/kinetics/", i_v, ".png"), plot = p1, width = 4, height = 3)
  }
}
#------------------------------------------------------------------------------#

#
# obtain MSD
#
{
  # computeMSD is the function from the library "flowcatchR."
  # I could not load it for some reasons, so pasted here as it is.
  computeMSD <- function(sx,sy,until=4){ 
    if(until > length(sx)-1){ print("Error: too long tau"); return(0)}
    msd.t <- rep(0,until)
    for (dt in 1:until){
      displacement.x <- as.vector(na.omit(sx[(1+dt):length(sx)]) - sx[1:(length(sx)-dt)])
      displacement.y <- as.vector(na.omit(sy[(1+dt):length(sy)]) - sy[1:(length(sy)-dt)])
      sqrdispl <- (displacement.x^2 + displacement.y^2)
      msd.t[dt] <- mean(sqrdispl)
    }
    return(msd.t)
  }
  
  load("data_fmt/df_expand_solo.rda")
  df_msd$msd <- NA
  for(i_v in unique(df_msd$video)){
    df_temp <- df_msd[df_msd$video == i_v,]
    msd <- computeMSD(df_temp$x, df_temp$y, dim(df_temp)[1]-1)
    df_msd[df_msd$video == i_v,]$msd <- c(NA, msd)
  }
  df_msd$time <- df_msd$frame/30
  dfMSD_solo <- df_msd
  
  load("data_fmt/df_expand_tandem.rda")
  
  df_msd$msd <- NA
  for(i_t in unique(df_msd$tandem_event)){
    #i_v <- unique(df_msd$video)[1]
    df_msd[df_msd$tandem_event == i_t, ]$frame <- df_msd[df_msd$tandem_event == i_t, ]$frame - df_msd[df_msd$tandem_event == i_t, ]$frame[1]
    df_msd[df_msd$tandem_event == i_t, ]$x <- df_msd[df_msd$tandem_event == i_t, ]$x - df_msd[df_msd$tandem_event == i_t, ]$x[1]
    df_msd[df_msd$tandem_event == i_t, ]$y <- df_msd[df_msd$tandem_event == i_t, ]$y - df_msd[df_msd$tandem_event == i_t, ]$y[1]
    df_temp <- subset(df_msd, tandem_event == i_t)
    if(df_temp$tandem[1] & dim(df_temp)[1]/5 > 200){
      msd <- computeMSD(df_temp$x, df_temp$y, dim(df_temp)[1]-1)
      df_msd[df_msd$tandem_event == i_t,]$msd <- c(NA, msd)
    }
  }
  df_msd <- na.omit(df_msd)
  df_msd$time <- df_msd$frame/30
  dfMSD_tandem <- df_msd
  
  load("data_fmt/df_antam.rda")
  df$msd <- NA
  for(i_t in unique(df$ind)){
    #i_t <- unique(df$ind)[1]
    df_temp <- subset(df, ind == i_t)
    msd <- computeMSD(df_temp$x, df_temp$y, dim(df_temp)[1]-1)
    df[df$ind == i_t,]$msd <- c(NA, msd)
  }
  dfMSD_ANTAM <- df
  
  save(dfMSD_tandem, dfMSD_solo, dfMSD_ANTAM, file = "data_fmt/df_msd.rda")
}
