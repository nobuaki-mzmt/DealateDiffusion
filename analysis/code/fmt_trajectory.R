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
  rad_near_wall <- dish_size/(2*sqrt(2))  # how far from the center of the wall to be near wall
  tandem_dis    <- 3.1 # ftip-mhead distance threshold based on Mizumoto and Reiter 2025 
  wall_following_thresh_dis <- 10 # mm
  angle_mismatch_threshold <- 0.5 # to investigate the mismatch between move_dir and head_dir

  circle_data <- data.frame(
    x = dish_size/2 * cos(seq(0, 2 * pi, length.out = 101)),
    y = dish_size/2 * sin(seq(0, 2 * pi, length.out = 101))
  )
  
  process_plot <- F
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# expand dish to open space
#------------------------------------------------------------------------------#
{
  ### Functions for edge corrections
  {
    unwrap_wall_following <- function(df_temp, process_plot = 0){
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
        
        if(duration > 1 & moved_dis > wall_following_thresh_dis){
          
          if(i_w == process_plot){
            p1 <- ggplot(data = df_temp[wall_frames,], aes(x = rad*cos(theta), y = rad*sin(theta)))+
              geom_path(color = "grey20", alpha = 0.5) +
              geom_point(data = df_temp[wall_frames[1],], aes(x = rad*cos(theta), y = rad*sin(theta)))+
              coord_fixed(xlim = c(-70, 70), ylim = c(-70, 70)) +
              scale_x_continuous(breaks = c(-70, 0, 70)) +
              scale_y_continuous(breaks = c(-70, 0, 70)) +
              theme_classic(base_family = "Arial") +
              theme(legend.position = "none") +
              labs(x = "X (mm)", y = "Y (mm)")
              #ggtitle(paste(i_w, "/", length(crockwise$lengths), round(sd(df_temp[wall_frames, ]$rad), 2)))
            ggsave(plot = p1, filename = paste0("output/", i_v, "_wraped_", i_w, ".pdf"), 
                   width = 3, height = 3, device = cairo_pdf)
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
            
            df_temp$head_dir[r_frames] <- df_temp$head_dir[r_frames] - rotation_per_frame[i]
            df_temp$move_dir[r_frames] <- df_temp$move_dir[r_frames] - rotation_per_frame[i]
          }
          
          if(i_w == process_plot){
            p2 <- ggplot(data = df_temp[wall_frames,], aes(x = Center_x, y = Center_y))+
              geom_path(color = "grey20", alpha = 0.5) +
              geom_point(data = df_temp[wall_frames[1],], aes(x = Center_x, y = Center_y))+
              coord_fixed() + 
              theme_classic(base_family = "Arial") +
              theme(legend.position = "none") +
              labs(x = "X (mm)", y = "Y (mm)")
            #ggtitle(paste(i_w, "/", length(crockwise$lengths), round(rotation))) #+ xlim(0,500) + ylim(0,500)
            ggsave(plot = p2, filename = paste0("output/", i_v, "_unwrap_", i_w, ".pdf"),
                   width = 3, height = 3, device = cairo_pdf)
          }
        } 
        step_count <- step_count + duration
      }
      return(df_temp)
    }
    
    # fix wall bounded turn
    # detect big turns near the edge and fix it
    correct_wall_bounded_turns <- function(df_temp){ 
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
      return(df_temp)
    }
    
    # look for it, based on the mismatch between head_dir and move_dir and fix it
    fix_directional_mismatch <- function(df_temp){
      df_temp$head_dir <- atan2(sin(df_temp$head_dir), cos(df_temp$head_dir))
      df_temp$move_dir <- atan2(sin(df_temp$move_dir), cos(df_temp$move_dir))
      angle_mismatch <- - df_temp$head_dir + df_temp$move_dir
      angle_mismatch <- atan2(sin(angle_mismatch), cos(angle_mismatch))
      
      mismatch_frame <- which(df_temp$rad > rad_near_wall)
      
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
        df_temp$move_dir[mismatch_frame[i]] <- df_temp$move_dir[mismatch_frame[i]] - rotation
        
        # add rest of the frame as the rotated place as an offset
        offset_pos <- as.numeric(df_temp1[mismatch_frame[i], 1:2])
        r_frames <- (mismatch_frame[i]+1):nrow(df_temp1)
        df_temp1[r_frames, 1] <- df_temp1[r_frames, 1] - offset_pos[1]
        df_temp1[r_frames, 2] <- df_temp1[r_frames, 2] - offset_pos[2]
        df_temp[r_frames, 1] <- df_temp1[r_frames, 1] + df_temp[mismatch_frame[i], 1]
        df_temp[r_frames, 2] <- df_temp1[r_frames, 2] + df_temp[mismatch_frame[i], 2]
      }
      
      df_temp$move_dir <- atan2(sin(df_temp$move_dir), cos(df_temp$move_dir))
      
      return(df_temp)
    }
    
  }
  
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
    
    df <- df %>% dplyr::select(Center_x, Center_y, rad, theta, head_dir, move_dir, step, video, colony, sex, frame)
    
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
      df_temp <- unwrap_wall_following(df_temp)
      
      # detect big turns near the edge and fix it
      #df_temp <- correct_wall_bounded_turns(df_temp)
      
      # fix wall bounded turn
      df_temp <- fix_directional_mismatch(df_temp)
      
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
    colnames(df)[colnames(df) %in% c("fCenter_x", "fCenter_y")] <- c("Center_x", "Center_y")

    df$rad <- sqrt(df$Center_x^2 + df$Center_y^2)
    df$theta <- atan2(df$Center_y, df$Center_x)
    
    df$head_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
    df$move_dir <- c(NA, atan2(diff(df$Center_y), diff(df$Center_x)))
    df$step <- c(NA, sqrt(diff(df$Center_y)^2 + diff(df$Center_x)^2))
    
    df <- df %>% dplyr::select(Center_x, Center_y, rad, theta, head_dir, move_dir, step, tandem, tandem_event, video, colony, frame)
    
    process_plot <- F
    
    for(i_v in unique(df$video)){
      
      print(i_v)
      df_temp <- df[df$video == i_v,]
      df_temp$step[1] <- NA
      
      p1 <- ggplot(df_temp, aes(x = Center_x, y = Center_y, col = tandem, group = tandem_event))+
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
        # convert circular movements to open area
        df_temp <- unwrap_wall_following(df_temp)
        
        # detect big turns near the edge and fix it
        #df_temp <- correct_wall_bounded_turns(df_temp)
        
        # fix wall bounded turn
        df_temp <- fix_directional_mismatch(df_temp)
      }
        
      # plot expanded trajectories
      {
        x_range <- range(df_temp$Center_x)
        y_range <- range(df_temp$Center_y)
        range_size <- max(diff(x_range), diff(y_range))
        center_xy <- c(mean(x_range), mean(y_range))
        x_limits <- c(center_xy[1] - range_size / 2 - dish_size, center_xy[1] + range_size / 2 + dish_size)
        y_limits <- c(center_xy[2] - range_size / 2 - dish_size, center_xy[2] + range_size / 2 + dish_size)
        
        p2 <- ggplot() +
          geom_path(data = df_temp, aes(x = Center_x, y = Center_y, 
                                        col = tandem, group = tandem_event),
                    alpha = 0.5, linewidth = 0.25) +
          geom_path(data = df_temp[!df_temp$tandem,], 
                    aes(x = Center_x, y = Center_y, group = tandem_event),
                    alpha = 1, linewidth = 0.4, col = viridis(2)[1]) +
          #geom_point(data = df_temp[!df_temp$tandem,], 
          #           aes(x = Center_x, y = Center_y),
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
      }
      
      df[df$video == i_v,] <- df_temp
    }
    
    
    df$acc  <- c(diff(df$step), NA)
    df$turn <- c(NA, diff(df$head_dir))
    df$turn <- atan2(sin(df$turn), cos(df$turn))
    df$turn[df$frame == 0] <- NA
    
    df_msd <- df %>% dplyr::select(frame, x = Center_x, y = Center_y, 
                                   tandem, tandem_event, video, colony)
    df <- df %>% dplyr::select(frame, head_dir = head_dir, move_dir = move_dir,
                               step = step, acc, turn, video, colony, tandem, tandem_event)
    
    save(df, df_msd, file = "data_fmt/df_expand_tandem.rda")
  }
  
  ### Output examples for figures
  {
    # data_prep
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
      colnames(df)[colnames(df) %in% c("fCenter_x", "fCenter_y")] <- c("Center_x", "Center_y")
      
      df$rad <- sqrt(df$Center_x^2 + df$Center_y^2)
      df$theta <- atan2(df$Center_y, df$Center_x)
      
      df$head_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
      df$move_dir <- c(NA, atan2(diff(df$Center_y), diff(df$Center_x)))
      df$step <- c(NA, sqrt(diff(df$Center_y)^2 + diff(df$Center_x)^2))
      
      df <- df %>% dplyr::select(Center_x, Center_y, rad, theta, head_dir, move_dir, step, tandem, tandem_event, video, colony, frame)
      
    }
    
    # unwrap
    {
      i_v = "Copfor_172_170627-01"
      # i_v = "Copfor_172_170627-05"
      
      df_temp <- df[df$video == i_v,]
      df_temp$step[1] <- NA
      
      # convert circular movements to open area
      df_temp <- unwrap_wall_following(df_temp, process_plot = 21)  
    }
    
    # discrepancy
    i_v = "Copfor_172_170627-10"
    df_temp <- df[df$video == i_v,]
    df_temp$step[1] <- NA
    df_temp$move_dir[1] <- NA
    
    p1 <- ggplot(data = df_temp[1:10,], aes(x = rad*cos(theta), y = rad*sin(theta)))+
      geom_path(color = "grey20", alpha = 0.5) +
      geom_point(data = df_temp[1,], aes(x = rad*cos(theta), y = rad*sin(theta)))+
      coord_fixed() +
      theme_classic() +
      theme(legend.position = "none") +
      labs(x = "X (mm)", y = "Y (mm)")
    p1
    ggsave(plot = p1, filename = paste0("output/", i_v, "_discrepancy_1-10", ".pdf"), 
           width = 3, height = 3, device = cairo_pdf)
    
    df_temp <- fix_directional_mismatch(df_temp)
    
    p2 <- ggplot(data = df_temp[1:10,], aes(x = Center_x, y = Center_y))+
      geom_path(color = "grey20", alpha = 0.5) +
      geom_point(data = df_temp[1,], aes(x = Center_x, y = Center_y))+
      coord_fixed() + 
      theme_classic() +
      theme(legend.position = "none") +
      labs(x = "X (mm)", y = "Y (mm)")
    ggsave(plot = p2, filename = paste0("output/", i_v, "_fixed_1-10", ".pdf"), 
           width = 3, height = 3, device = cairo_pdf)
    
    
    
    
  }
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
  
  save(dfMSD_tandem, dfMSD_solo, file = "data_fmt/df_msd.rda")
}
