# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 22:11:59 2024

@author: ramne
"""

from collections import defaultdict
import numpy as np
from skimage.io import imshow
from skimage.measure import label
from scipy.ndimage.morphology import distance_transform_edt
import matplotlib.pyplot as plt
from PIL import Image
import PIL
import glob

import imageio.v3 as iio
import skimage as ski


def unet_weight_map(y, wc=None, w0 = 10, sigma =5):

    """
    Generate weight maps as specified in the U-Net paper
    for boolean mask.
    
    "U-Net: Convolutional Networks for Biomedical Image Segmentation"
    https://arxiv.org/pdf/1505.04597.pdf
    
    Parameters
    ----------
    mask: Numpy array
        2D array of shape (image_height, image_width) representing binary mask
        of objects.
    wc: dict
        Dictionary of weight classes.
    w0: int
        Border weight parameter.
    sigma: int
        Border width parameter.
    Returns
    -------
    Numpy array
        Training weights. A 2D array of shape (image_height, image_width).
    """
    y = y.reshape(y.shape[0], y.shape[1])
    print(y)
    labels = label(y)
    print(labels)
    no_labels = labels == 0
    print(no_labels)
    label_ids = sorted(np.unique(labels))[1:]
    print(label_ids)
    

    if len(label_ids) > 1:
        distances = np.zeros((y.shape[0], y.shape[1], len(label_ids)))

        for i, label_id in enumerate(label_ids):
            distances[:,:,i] = distance_transform_edt(labels != label_id)

        distances = np.sort(distances, axis=2)
        d1 = distances[:,:,0]
        d2 = distances[:,:,1]
        w = w0 * np.exp(-1/2*((d1 + d2) / sigma)**2) * no_labels
        
        if wc:
            class_weights = np.zeros_like(y)
            for k, v in wc.items():
                class_weights[y == k] = v
            w = w + class_weights
    else:
        w = np.zeros_like(y)
    
    return w 

# upload de binary mask that you have made.
y = Image.open(r"C:\Users\ramne\OneDrive\Documenten\Master BFW jaar 1\Reseach project 1\40x images muriel zonder rgba channel\binary_watershed\input01_40.tif")
y= np.array(y)
print(y)


wc = {
    0: 0, # background
    1: 1  # objectsplt.show()
}

w = unet_weight_map(y, wc)

plt.imshow(w,'gray')
plt.axis('off')
plt.savefig('wmap.tif',bbox_inches ="tight",pad_inches=0)# name it how you want it
plt.show()

#resize wmap to 512x512 pixels.
width = 256
height = 256
image=Image.open(r"C:\Users\ramne\OneDrive\Documenten\Master BFW jaar 1\Reseach project 1\Deepsea\weight map python\wmap.tif") #change name of how you save the wmap plot.
image=image.resize((width,height), PIL.Image.LANCZOS)
image.save('resized_wmap.tif')
print(image.size) #512x512 

shapes01 = iio.imread(uri=r"C:\Users\ramne\OneDrive\Documenten\Master BFW jaar 1\Reseach project 1\Deepsea\weight map python\resized_wmap.tif")[:,:,:3]
fig, ax = plt.subplots()
plt.axis('off')
plt.imshow(shapes01)
# convert the image to grayscale

gray_shapes = ski.color.rgb2gray(shapes01)

# blur the image to denoise
blurred_shapes = ski.filters.gaussian(gray_shapes, sigma=1.0)

fig, ax = plt.subplots()
plt.imshow(blurred_shapes, cmap="gray")
plt.axis('off')

# create a histogram of the blurred grayscale image
histogram, bin_edges = np.histogram(blurred_shapes, bins=256, range=(0.0, 1.0))

fig, ax = plt.subplots()
plt.plot(bin_edges[0:-1], histogram)
plt.title("Grayscale Histogram")
plt.xlabel("grayscale value")
plt.ylabel("pixels")
plt.xlim(0, 1.0)

# create a mask based on the threshold
t = 0.005
binary_mask = blurred_shapes > t

fig, ax = plt.subplots()
plt.imshow(binary_mask, cmap="gray")
plt.axis('off')#
plt.savefig('binary_wmap.tif',bbox_inches ="tight",pad_inches=0)# name it how you want it
plt.show()

#resize binary wmap
width = 256
height = 256
image=Image.open(r"C:\Users\ramne\OneDrive\Documenten\Master BFW jaar 1\Reseach project 1\Deepsea\weight map python\binary_wmap.tif")
image=image.resize((width,height), PIL.Image.LANCZOS)
image.save('resized_binary_wmap.tif')
print(image.size) #512x512 


