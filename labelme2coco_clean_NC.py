# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 20:16:45 2019

@author: jizh
"""
# adapted from: 
# http://www.programmersought.com/article/22821169839/;jsessionid=A385E503DE5519BA214D561ABD574B10

###                      prepare the csv file 

## use labelme to label images manually

###############################################################################
import argparse
import json
import matplotlib.pyplot as plt
import skimage.io as io
import cv2
from labelme import utils
import numpy as np
import glob
import PIL.Image

#
def image(data,num):
        image={}
        img = utils.img_b64_to_arr(data['imageData'])  
        height, width = img.shape[:2]
        img = None
        image['height']=height
        image['width'] = width
        image['id']=num+1
        image['file_name'] = data['imagePath'].split('/')[-1]
        height=height
        width=width
        return image

def annotation(points,label,num):
        annotation={}
        annotation['iscrowd'] = 0
        annotation['image_id'] = num+1

        annotation['bbox'] = list(map(float, getbbox(points)))

        annotation['category_id'] = getcatid(label)
        annotation['id'] = annID
        return annotation

def getbbox(points):
        polygons = points
        mask = polygons_to_mask([height,width], polygons)
        return mask2box(mask)

def getcatid(label):
        for categorie in categories:
            if label[0]==categorie['name']:
                return categorie['id']
        return -1


def polygons_to_mask(img_shape, polygons):
        mask = np.zeros(img_shape, dtype=np.uint8)
        mask = PIL.Image.fromarray(mask)
        xy = list(map(tuple, polygons))
        PIL.ImageDraw.Draw(mask).polygon(xy=xy, outline=1, fill=1)
        mask = np.array(mask, dtype=bool)
        return mask


def mask2box( mask):
        index = np.argwhere(mask == 1)
        rows = index[:, 0]
        clos = index[:, 1]
        # 
        left_top_r = np.min(rows)
        left_top_c = np.min(clos)#.tolist()  # x

        right_bottom_r = np.max(rows) #.tolist()
        right_bottom_c = np.max(clos)#.tolist()

        return [left_top_c, left_top_r, right_bottom_c, right_bottom_r]  # 


def categorie(label):
        categorie={}
        categorie['supercategory'] = label[0]
        categorie['id']=len(label)+1 # 
        categorie['name'] = label[0]
        return categorie

def concatenate_list_data(x):
    result= list()
    for element in x:
        result.append(str(int(element)))
    return result

#############
# extract information from json file
from tqdm import tqdm
import os

path = "path to folder of json files"
train_ids = next(os.walk(path))[2]
train_ids

# initialize all arguments
images=[]
categories=[]
annotations=[]
label=[]
height=0
width=0

for n, id_ in tqdm(enumerate(train_ids), total=len(train_ids)):
    print(id_)
    
    num=int(id_.split(".")[0])
    json_file= os.path.join(path,id_)
    with open(json_file,'r') as fp:
        data = json.load(fp)
        annID=1
        
        image_x=image(data,num)
        images.append(image_x)
        height=image_x["height"]
        width=image_x["width"]
        for shapes in data['shapes']:
            label=shapes['label'].split('_')
            print(label)
            if label[0] not in label:
                categories.append(categorie(label))
                label.append(label[0])
            points=shapes['points']
            annotations.append(annotation(points,label,num))
            annID+=1


data_tomato={}
data_tomato['images']=images
data_tomato['categories']=categories
data_tomato['annotations']=annotations

#################################################################################

#########
labels={"1": "tomato"}
ann_img_list=list()
for i in range(len(data_tomato["images"])):
    print("i", i)
    img_id_i=data_tomato["images"][i]["file_name"]
    img_id=int(img_id_i.split(".")[0])
    for x in range(len(data_tomato["annotations"])):
            if data_tomato["annotations"][x]["image_id"]-1 ==img_id:
                bbox_i=data_tomato["annotations"][x]["bbox"]
                lbel_i=labels["1"]
                bbox_i_l=concatenate_list_data(bbox_i)
                ann_img_i=",".join(bbox_i_l)+","+lbel_i
                ann_img=img_id_i+"," +ann_img_i
                
                x=x+1
#               print(ann_img)
                ann_img_list.append(ann_img)

with open('annotation_tomato.txt', 'w') as f:
    for item in ann_img_list:
        f.write("%s\n" % item)






