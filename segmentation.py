import subprocess
import os
import numpy as np
import cv2
from skimage.measure import label, regionprops
import matplotlib.pyplot as plt
import pandas as pd
import imageio
import imagej
import scipy.ndimage
import skimage.filters
import xlsxwriter
import sys
import openpyxl
import scyjava as sj
import csv
import cv2
import numpy as np
from scipy.spatial import distance
import numpy as np
from skimage import io, filters
import tifffile as tiff
from scipy.ndimage import binary_fill_holes
import numpy as np
import skimage as ski
import scipy as sp
from skimage.exposure import histogram
import os
import matplotlib.pyplot as plt
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, threshold_niblack, threshold_sauvola
from matplotlib.colors import ListedColormap
from skimage.util import montage as skimage_montage
from scipy.ndimage import distance_transform_edt
from skimage import color, morphology
from scipy.ndimage import label, sum as ndi_sum

file = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\Set 1\mirco_gell.tif"
test_file = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\Set 1\MAX_mirco_gell.tif"
segmented_image = r"C:\Users\aryan\PycharmProjects\pythonProject\.venv\Scripts\updated_segmented_multislice_colored_col.tiff"
image_3d = tiff.imread(file)
image_2d = tiff.imread(test_file)

def filter_otsu_3d(image):
    otsu_thres = filters.threshold_otsu(image)
    bimage = np.zeros_like(image)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            for k in range(image.shape[2]):
                if image[i, j, k] > otsu_thres:
                    bimage[i, j, k] = 1
                else:
                    bimage[i, j, k] = 0
    return bimage

def filter_otsu_2d(image):
    otsu_thres = filters.threshold_otsu(image)
    bimage = np.zeros_like(image)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            if image[i, j] > otsu_thres:
                bimage[i, j] = 1
            else:
                bimage[i, j] = 0
    return bimage


binary_image = filter_otsu_3d(image_3d)

distance = ndi.distance_transform_edt(binary_image)

local_max_coords = ski.feature.peak_local_max(distance, min_distance=20, exclude_border=True)
local_max_mask = np.zeros(distance.shape, dtype=bool)
local_max_mask[tuple(local_max_coords.T)] = True
markers = ski.measure.label(local_max_mask)
colors = np.array([
    [255, 0, 0]
])
segmented_cells = ski.segmentation.watershed(-distance, markers, mask=binary_image)
segmented_cells = ski.color.label2rgb(segmented_cells,colors=colors, bg_label=0)


labeled_array, num_features = label(segmented_cells)

# Calculate volumes
volumes = ndi_sum(segmented_cells > 0, labeled_array, range(1, num_features + 1), )
volumes = volumes * 1.089451638
print("Labeled Array:")
print(labeled_array)
print("Number of Features:", num_features)
print("Volumes of Regions:", sorted(volumes))

for i in sorted(volumes):
    print(i)

output_tiff_file = 'updated_segmented_multislice_colored.tiff'
with tiff.TiffWriter(output_tiff_file, bigtiff=True) as tiff_writer:
    for i in range(segmented_cells.shape[0]):
        tiff_writer.write(segmented_cells[i])

print(f"Segmented images saved to {output_tiff_file}")