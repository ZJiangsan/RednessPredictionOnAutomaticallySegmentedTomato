The scripts here aim to predict redness of automatically segmented tomato through hyprespectral image reconstruction from single RGB image

# keras_resnet_tomato_detection_segmentation
Google Colab version adapted from  https://github.com/fizyr/keras-retinanet
can go follow the link directly if there is anything unclear
# HSCNNR_tomato_rgb2HSI_NC_TSsplit_noAugm
python code can be run on Google Colab. The major architecture are adapted from HSCNNR, http://openaccess.thecvf.com/content_cvpr_2018_workshops/papers/w13/Shi_HSCNN_Advanced_CNN-Based_CVPR_2018_paper.pdf.
Data processing and model evaluation metrics are adapted from https://github.com/SaoYan/Multiscale-Super-Spectral
# labelme2coco
This file converts json files generated from labelme to csv files of image names and bounding boxes of tomatoes inside image
# annocation_tomato_segmentation_NC
The annotations of these tomatoes in the images: first column lists original images's file path while the last column includes the path of these masks of the original images; four columns in the middle indicate where the bounding boxes of tomatoes are.
# class_tomato
Show one example class name in csv file
