data_name <- paste0("data_fmt/tandem_df.feather")
df <- arrow::read_feather(data_name)
df <- as.data.frame(df)

df[,3:14] <- df[,3:14] - 70

df <- df[df$video == "Copfor_172_170627-01",]

f_dir <- atan2(df$fHead_y - df$fTip_y, df$fHead_x - df$fTip_x)
f_move_dir <- atan2(c(diff(df$fCenter_y), NA), c(diff(df$fCenter_x), NA))
f_move_dis <- c(sqrt(diff(df$fCenter_y)^2 + diff(df$fCenter_x)^2), NA)

fx <- df$fCenter_x
fy <- df$fCenter_y

for(i in 2:length(fx)){
  fx[i] <- fx[i-1] + f_move_dis[i-1] * cos(f_dir[i-1])
  fy[i] <- fy[i-1] + f_move_dis[i-1] * sin(f_dir[i-1])
}


plot(fx[830:890], fy[830:890], type ="l")
plot(df$fCenter_x[830:890], df$fCenter_y[830:890], type ="l")


df$fCenter_x <- fx
df$fCenter_y <- fy

df$f_r <- sqrt(fx^2 + fy^2)
df$wall <- df$f_r > 30#70/sqrt(2)
df$theta <- atan2(fy, fx)


df_temp <- df[df$video == "Copfor_172_170627-01", c(11:12, 18:20)]

ggplot(df_temp, aes(x = fCenter_x, y = fCenter_y))+
  #geom_point(alpha = 0.1) + 
  geom_path(alpha = 0.5) +
  xlim(-70,70)+ylim(-70,70)+theme(aspect.ratio = 1)

crockwise <- c(0, diff(df_temp$theta)) < 0
crockwise <- rle(crockwise)



step_count <- 0
for(i_w in 1:length(crockwise$lengths)){
  duration <- crockwise$length[i_w]
  
  #if(wall_follow$values[i_w]){
  if(duration > 5){
    
    wall_frames <- 1:duration + step_count
    
    p1 <- ggplot(data = df_temp[wall_frames,], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
      geom_path(color = 1, alpha = 0.5) +
      geom_point(data = df_temp[wall_frames[1],], aes(x = f_r*cos(theta), y = f_r*sin(theta)))+
      xlim(-70,70) + ylim(-70,70) + 
      coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(sd(df_temp[wall_frames, ]$f_r), 2))) #+ xlim(0,500) + ylim(0,500)
    
    rotation <- sum( diff(df_temp[wall_frames, ]$theta) )
    rotation_per_frame <- rotation/length(diff(df_temp[wall_frames, ]$theta))
    
    for(i in wall_frames[-1]){
      df_temp2 <- df_temp
      
      df_temp2[i:nrow(df_temp2), 1] <- df_temp2[i:nrow(df_temp2), 1] - df_temp[i-1, 1]
      df_temp2[i:nrow(df_temp2), 2] <- df_temp2[i:nrow(df_temp2), 2] - df_temp[i-1, 2]
      
      df_temp[i:nrow(df_temp), 1] <- cos(rotation_per_frame)*df_temp2[i:nrow(df_temp2), 1] + sin(rotation_per_frame)*df_temp2[i:nrow(df_temp2), 2] 
      df_temp[i:nrow(df_temp), 2] <- cos(rotation_per_frame)*df_temp2[i:nrow(df_temp2), 2] - sin(rotation_per_frame)*df_temp2[i:nrow(df_temp2), 1] 
      
      df_temp[i:nrow(df_temp), 1] <- df_temp[i:nrow(df_temp), 1] + df_temp[i-1, 1]
      df_temp[i:nrow(df_temp), 2] <- df_temp[i:nrow(df_temp), 2] + df_temp[i-1, 2]
    }
    
    p2 <- ggplot(data = df_temp[wall_frames,], aes(x = fCenter_x, y = fCenter_y))+
      geom_path(color = 1, alpha = 0.5) +
      geom_point(data = df_temp[wall_frames[1],], aes(x = fCenter_x, y = fCenter_y))+
      coord_fixed() + ggtitle(paste(i_w, "/", length(crockwise$lengths), round(rotation))) #+ xlim(0,500) + ylim(0,500)
    
    #ggsave(plot = p1+p2, filename = paste0(i_w, ".png"))
    
  } 
  
  step_count <- step_count + duration
  
  if(i_w %% 100 == 101){
    p <- ggplot()+
      geom_path(data = df_temp[(step_count+1):nrow(df_temp),], aes(x = fCenter_x, y = fCenter_y), color = 1, alpha = 0.5) +
      geom_path(data = df_temp[1:step_count,], aes(x = fCenter_x, y = fCenter_y), color = 2, alpha = 0.5) +
      coord_fixed() + ggtitle(paste(i_w, "/", length(wall_follow$lengths))) #+ xlim(0,500) + ylim(0,500)
    print(p)
    Sys.sleep(0.1)
  }
}

ggplot(df_temp, aes(x = fCenter_x, y = fCenter_y))+
  #geom_point(alpha = 0.1) + 
  geom_path(alpha = 0.5) +
  coord_fixed()
