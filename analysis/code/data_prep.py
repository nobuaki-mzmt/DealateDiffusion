"""
data_preparation.py
N. Mizumoto
This script reads all .h5 results from SLEAP and organize for the further analysis

Cop_for
"""
import sys
import os
import glob

import h5py

import numpy as np
from numpy.linalg import norm

import pandas as pd

import scipy
from scipy.interpolate import interp1d

import math
import feather

#------------------------------------------------------------------------------#
# interpolate the data
#------------------------------------------------------------------------------#
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
        if len(x) > 3:
          f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
          # Fill missing
          xq = np.flatnonzero(np.isnan(y))
          y[xq] = f(xq)
          # Fill leading or trailing NaNs with the nearest non-NaN values
          mask = np.isnan(y)
          y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
          Y[:, i] = y
          if sum(np.isnan(y)) > 0:
            print("error"+str(i))
            print("error"+(i))
    # Restore to initial shape.
    Y = Y.reshape(initial_shape)
    return Y
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main():
  dish_size = 145
  window = 5
  weights = np.ones(window) / window
  pad_width = (window - 1) // 2

  ### 1. tandem data processing
  if(False):
    print("tandem data processing")
    data_place = glob.glob("data_raw/tandem/*")
    
    df = pd.DataFrame()
    skip_list = []  # this is for debug
    
    for f_name in data_place:
      
      print(f_name)
      
      ## metadata
      pair_name = os.path.basename(f_name.replace(".h5", ""))
      species = pair_name.split("_")[0]
      colony = pair_name.split("_")[1]
      
      ## load data
      with h5py.File(f_name, "r") as f:
        try:
          locations = f['tracks'][:]
          node_names = [n.decode() for n in f["node_names"][:]]
          print("'tracks' exists")
        except KeyError:
          print("'tracks' does not exist")
          skip_list.append(f_name)
          continue;
      
      locations = locations.T
      total_frame = locations.shape[0]
  
      # remove tracks with > 2
      if locations.shape[3] > 2:
        print("there are too many tracks")
        locations = locations[:,:,:,0:2]
    
      ## processing locations
      # data filling
      locations = fill_missing(locations)
      
      # filtering
      for i_ind in range(locations.shape[3]):
        for i_coord in range(locations.shape[2]):
          for i_nodes in range(locations.shape[1]):
            locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
            locations[:, i_nodes, i_coord, i_ind] = np.convolve( np.pad(locations[:, i_nodes, i_coord, i_ind], pad_width, mode='edge'), weights, mode='valid')
            
      # scaling in mm (1200 pixels = dish_size)
      locations[:, :, :, :] = locations[:, :, :, :] / 1200 * dish_size
  
      # summarize data
      df_temp = {
          "frame": list(range(0,locations.shape[0],1)),
          "fHead_x": locations[:, node_names.index('headtip'), 0, 0].round(2),
          "fHead_y": locations[:, node_names.index('headtip'), 1, 0].round(2),
          "mHead_x": locations[:, node_names.index('headtip'), 0, 1].round(2),
          "mHead_y": locations[:, node_names.index('headtip'), 1, 1].round(2),
          "fTip_x": locations[:, node_names.index('abdomentip'), 0, 0].round(2),
          "fTip_y": locations[:, node_names.index('abdomentip'), 1, 0].round(2),
          "mTip_x": locations[:, node_names.index('abdomentip'), 0, 1].round(2),
          "mTip_y": locations[:, node_names.index('abdomentip'), 1, 1].round(2),
          #"fCenter_x": locations[:, node_names.index('abdomenfront'), 0, 0].round(2),
          #"fCenter_y": locations[:, node_names.index('abdomenfront'), 1, 0].round(2),
          #"mCenter_x": locations[:, node_names.index('abdomenfront'), 0, 1].round(2),
          #"mCenter_y": locations[:, node_names.index('abdomenfront'), 1, 1].round(2),
          "fCenter_x": locations[:, node_names.index('headtip'), 0, 0].round(2) + (locations[:, node_names.index('headtip'), 0, 0].round(2) - locations[:, node_names.index('abdomentip'), 0, 0].round(2))/2,
          "fCenter_y": locations[:, node_names.index('headtip'), 1, 0].round(2) + (locations[:, node_names.index('headtip'), 1, 0].round(2) - locations[:, node_names.index('abdomentip'), 1, 0].round(2))/2,
          "mCenter_x": locations[:, node_names.index('headtip'), 0, 1].round(2) + (locations[:, node_names.index('headtip'), 0, 1].round(2) - locations[:, node_names.index('abdomentip'), 0, 1].round(2))/2,
          "mCenter_y": locations[:, node_names.index('headtip'), 1, 1].round(2) + (locations[:, node_names.index('headtip'), 1, 1].round(2) - locations[:, node_names.index('abdomentip'), 1, 1].round(2))/2,
          
          "video": pair_name,
          "species": species,
          "colony": species + colony
        }
      df_temp = pd.DataFrame(df_temp)
      df_temp['group'] = df_temp['frame'] // 6
      downsampled_df_temp = df_temp.groupby('group').agg({
        'frame': 'first',
        'fHead_x': 'mean',
        'fHead_y': 'mean',
        'mHead_x': 'mean',
        'mHead_y': 'mean',
        'fTip_x': 'mean',
        'fTip_y': 'mean',
        'mTip_x': 'mean',
        'mTip_y': 'mean',
        'fCenter_x': 'mean',
        'fCenter_y': 'mean',
        'mCenter_x': 'mean',
        'mCenter_y': 'mean',
        "video": 'first',
        "species": 'first',
        "colony": 'first'
        }).reset_index(drop=True)
      df = pd.concat([df, pd.DataFrame(downsampled_df_temp)])
    
    #df = df[df["frame"] % 6 == 0]
    
    df.reset_index().to_feather("data_fmt/tandem_df.feather")
    print("skipped item: " + str(skip_list))
    
  ### 2. solo data processing
  if(True):
    print("solo data processing")
    data_place = glob.glob("data_raw/solo/*")
    df_dish = pd.read_csv("data_raw/df_dishsize_solo.csv")
  
    df = pd.DataFrame()
    skip_list = []  # this is for debug
    
    for f_name in data_place:
      
      print(f_name)
      
      ## metadata
      video_name = os.path.basename(f_name.replace(".h5", ""))
      colony = video_name.split("_")[1]
      sex    = video_name.split("_")[2]
      
      ## load data
      with h5py.File(f_name, "r") as f:
        try:
          locations = f['tracks'][:]
          node_names = [n.decode() for n in f["node_names"][:]]
          print("'tracks' exists")
        except KeyError:
          print("'tracks' does not exist")
          skip_list.append(f_name)
          continue;
      
      locations = locations.T
      total_frame = locations.shape[0]
  
      # remove tracks with > 1
      if locations.shape[3] > 1:
        print("there are too many tracks")
        locations = locations[:,:,:,0:1]
    
      ## processing locations
      # data filling
      locations = fill_missing(locations)
      
      # filtering
      for i_ind in range(locations.shape[3]):
        for i_coord in range(locations.shape[2]):
          for i_nodes in range(locations.shape[1]):
            locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
            #print(locations[0:10, i_nodes, i_coord, i_ind])
            locations[:, i_nodes, i_coord, i_ind] = np.convolve( np.pad(locations[:, i_nodes, i_coord, i_ind], pad_width, mode='edge'), weights, mode='valid')
            #print(locations[0:10, i_nodes, i_coord, i_ind])
      
      # scaling in mm
      dish_area = df_dish[df_dish['video'] == video_name].drop(columns=['video']).values.flatten()
      locations[:, :, 0, :] = (locations[:, :, 0, :] - dish_area[0]) / (dish_area[2]-dish_area[0]) * dish_size - dish_size/2
      locations[:, :, 1, :] = (locations[:, :, 1, :] - dish_area[1]) / (dish_area[3]-dish_area[1]) * dish_size - dish_size/2
  
      df_temp = {
          "frame": list(range(0,locations.shape[0],1)),
          "Head_x": locations[:, node_names.index('headtip'), 0, 0].round(2),
          "Head_y": locations[:, node_names.index('headtip'), 1, 0].round(2),
          "Tip_x": locations[:, node_names.index('abdomentip'), 0, 0].round(2),
          "Tip_y": locations[:, node_names.index('abdomentip'), 1, 0].round(2),
          #"Center_x": locations[:, node_names.index('abdomenfront'), 0, 0].round(2),
          #"Center_y": locations[:, node_names.index('abdomenfront'), 1, 0].round(2),
          "Center_x": locations[:, node_names.index('headtip'), 0, 0].round(2) + (locations[:, node_names.index('headtip'), 0, 0].round(2) - locations[:, node_names.index('abdomentip'), 0, 0].round(2))/2,
          "Center_y": locations[:, node_names.index('headtip'), 1, 0].round(2) + (locations[:, node_names.index('headtip'), 1, 0].round(2) - locations[:, node_names.index('abdomentip'), 1, 0].round(2))/2,
          "video": video_name,
          "colony": "Copfor" + colony,
          "sex": sex
          }
      df_temp = pd.DataFrame(df_temp)
      df_temp['group'] = df_temp['frame'] // 6
      downsampled_df_temp = df_temp.groupby('group').agg({
        'frame': 'first',
        'Head_x': 'mean',
        'Head_y': 'mean',
        'Tip_x': 'mean',
        'Tip_y': 'mean',
        'Center_x': 'mean',
        'Center_y': 'mean',
        "video": 'first',
        "colony": 'first',
        "sex": 'first'
        }).reset_index(drop=True)
      
      
      df = pd.concat([df, pd.DataFrame(downsampled_df_temp)])
    
    #df = df[df["frame"] % 6 == 0] # down sample to 5FPS (4.995 FPS)
    
    df.reset_index().to_feather("data_fmt/solo_df.feather")
    print("skipped item: " + str(skip_list))
    
  return 0

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
  # Call the main function without any parameter
  main()
#------------------------------------------------------------------------------#



