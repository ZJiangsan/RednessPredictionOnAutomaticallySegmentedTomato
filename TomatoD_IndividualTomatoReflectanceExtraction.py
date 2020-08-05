#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 23:00:33 2020

@author: zhaojiangsan
"""
##     TomatoD_HSI_RGB_extraction

import os
import h5py
import cv2
import numpy as np

from sklearn.model_selection import train_test_split
import hdf5storage

def save_matv73(mat_name, var_name, var):
    hdf5storage.savemat(mat_name, {var_name: var}, format='7.3', store_python_metadata=True)
    
path_hsi_fd = "/Users/zhaojiangsan/Documents/RScript/tomato_hsi_reconstructed/" # folder of hsi files

path_mask_fd = path_hsi_fd.replace('hsi_reconstructed', 'mask') # folder of mask files
var_ids = next(os.walk(path_mask_fd))[1] # 


for var_id in var_ids[1:]: 
    print("var_id", var_id)
    path_rgb = os.path.join(path_mask_fd, var_id)
    rgb_ids = next(os.walk(path_rgb))[2]

    for rgb_id in rgb_ids: # three rgb images for variety D
        print("rgb_id", rgb_id)
        rgb_id_path = os.path.join(path_rgb, rgb_id)
        ck_id = "_".join(((rgb_id.split(".")[0]).split("/")[-1]).split("_")[-3:]) # the original file name
        print("ck_id", ck_id)
        rgb_redR =  cv2.imread(rgb_id_path) # read the rgb file
        
        hsi_var__i = np.empty((0,205)) # create an empty np array 
        # create 4 patches for corresponding to reconstructed HSIs
        for m in range(2):
                for n in range(2):
                    rgb_x = rgb_redR[378*4*m:378*4*(m+1), 504*4*n:504*4*(n+1)]
                    
                    ## read reconstructed hsi
                    path_hsi = os.path.join(path_hsi_fd, var_id)
                    hsi_id = "{}_{}_{}.mat".format(ck_id, m, n) # hsi file to be read
                    
                    print("hsi_id", hsi_id)
                    
                    mat_i =  h5py.File(os.path.join(path_hsi, hsi_id),'r')
                    hsi_i =  np.float32(np.array(mat_i['img']))
                    # alternate the dimension to be similar like the rgb image (H*W*C)
                    hsi_j = np.transpose(hsi_i, [2,1,0])
                    
                    mat_i.close()
                        
                    
    
                    for p in list(np.unique(rgb_x[:,:, 0]))[1:]: # here to hsi of each single tomato
                        print("tomato_id:", int((p-16)/20))
                        
                        ## directly subsample the data and save due to the large size data
                        shit_x =  hsi_j[rgb_x[:,:,0]==int(p)]
                        shit_y = (np.ones(shit_x.shape[0])*((p-16)/20)).astype(int)
                        shit_x_save, shit_x_test, shit_y_save, shit_y_test = train_test_split(shit_x, shit_y, test_size = 0.7, random_state = 0)
                        
                        shit_x_y_save = np.concatenate((shit_y_save[:, None], shit_x_save), axis = 1)
                        hsi_var__i = np.concatenate((hsi_var__i, shit_x_y_save), axis = 0)
                        print("hsi_var__i", hsi_var__i.shape)
        # save the data in .mat file
        var_name = 'img' # the variable name used in mat file
        
        mat_name = ck_id + '.mat'
        mat_path = path_rgb.replace("_mask", "_hsi_mat")
        if not os.path.exists(mat_path):
            os.makedirs(mat_path)
        mat_save_path = os.path.join(mat_path,mat_name)
        save_matv73(mat_save_path, var_name, hsi_var__i)






