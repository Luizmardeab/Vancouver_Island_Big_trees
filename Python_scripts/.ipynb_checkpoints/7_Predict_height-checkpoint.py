import rasterio
from rasterio.transform import from_origin
import os
import numpy as np
import random
import glob
import tifffile as tiff
import os
import numpy as np
import sys
import tensorflow as tf
from scipy.ndimage import rotate
from PIL import Image
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.optimizers.schedules import CosineDecay
from tensorflow.keras.callbacks import CSVLogger, ModelCheckpoint, EarlyStopping
from tensorflow.keras.preprocessing.image import ImageDataGenerator
sys.path.append('/home/barros/Python_scripts/')
from smooth_tile_predition import predict_img_with_smooth_windowing
from Deep_learning_UNET import get_unet_128
from matplotlib import pyplot as plt
# Define function to extract unique identifier from the file name
import re

def extract_identifier(filename):
    match = re.search(r'\d+_\d+_\d+', filename)
    if match:
        return match.group(0)
    return None
import gc

root_directory = '/home/barros/Big_tree/'
os.chdir(root_directory)

image_dir = os.path.join(root_directory, "sentile4predVSG/")
mask_dir = os.path.join(root_directory, "Masks/")

numvar = ['009','010', '011', '012', '020', '021', '022', '023', '024', 
           '031', '032', '033', '034', '035', '036', '042', '043', 
          '044', '045', '046', '047','048', '054', '055', '056', 
          '057', '058', '059', '060', '065', '066', '067', '068', 
          '069', '070', '071', '072', '076', '077', '078', '079', 
          '080', '081', '082', '083', '084', '087', '088', '089', 
          '090', '091', '092', '093', '094', '095', '096', '099', 
           '100', '101', '102', '103', '104', '105', '106', '107', 
          '108', '109', '110', '111', '112', '113', '114', '115', 
           '116', '117', '118', '119', '120', '121', '122', '123', 
          '124', '125', '126', '127', '128', '129', '130', '131', '132']

          
# Define the pattern to match tile files (e.g., 'tile_*.tif')
tile_pattern = os.path.join(image_dir, 'Tile_*.tif') #Tile_ changed from Tile_

# Use glob to find all tile file paths matching the pattern
tile_paths = glob.glob(tile_pattern)

# Filter paths to include only the desired numbers
filtered_image_paths = [tile_path for tile_path in tile_paths if any(num in tile_path for num in numvar)]

filtered_image_paths = sorted(filtered_image_paths)
# Print the list of filtered tile file paths
print("Filtered tile file paths:")
for tile_path in filtered_image_paths:
    print(tile_path)

###########################################################################

IMG_HEIGHT = 256#image.shape[1]
IMG_WIDTH  = 256#image.shape[2]
IMG_CHANNELS = 9#image.shape[3]+1
n_classes=1

print("height= {}, Width = {}, Channels= {}, and classes = {}".format(IMG_HEIGHT,IMG_WIDTH,IMG_CHANNELS,n_classes))


output_pred_dir = os.path.join(root_directory, "pred_S12_10mGLOBRMSE/")

# Create the directory if it doesn't exist
os.makedirs(output_pred_dir, exist_ok=True)

# size of patches
patch_size = 256

# Number of classes 
n_classes = 1

# List of best model weight files
best_model_paths = [
    'Model/Sent12_10_GlobalHisErrMAE/unet_tf2_143_1.0892.h5',
    'Model/Sent12_10_GlobalHisErrMAE/unet_tf2_139_1.1094.h5',
    'Model/Sent12_10_GlobalHisErrMAE/unet_tf2_122_1.1304.h5',
    'Model/Sent12_10_GlobalHisErrMAE/unet_tf2_114_1.1437.h5',
    'Model/Sent12_10_GlobalHisErrMAE/unet_tf2_098_0.0688.h5'
]

# Compile the model
def get_model():
    return get_unet_128(n_classes=n_classes, IMG_HEIGHT=IMG_HEIGHT, IMG_WIDTH=IMG_WIDTH, IMG_CHANNELS=IMG_CHANNELS, activation = activation_function)

def normalize_band(band, mean_val, std_val):
    # Normalize the band so that mean is 0, max is 1, and min is -1
    norm_band = (band - mean_val) / std_val  # Standardize
    return np.clip(norm_band, -6, 6)  # Clip to ensure the values stay between -1 and 1
    
for image_path in filtered_image_paths[20:]:
    # Initialize array to accumulate predictions
    all_preds = []
    # Predict with each model and accumulate
    for model_path in best_model_paths:
        model = get_model()
        model.compile(optimizer=optimizer,
                    loss=laplace_loss, 
                metrics= metrics)
        # model.summary()
        model.load_weights(os.path.join(root_directory,model_path))

        image = tiff.imread(image_path)
        # image = np.transpose(image, (1, 2, 0))
        with rasterio.open(image_path) as src:
            original_transform = src.transform
            original_crs = src.crs
            original_meta = src.meta

        # Extract the number using regular expressions (assuming the number is two digits)
        number = re.search(r'_(\d{3})\.', image_path).group(1)
        print(f"Opened image: {number}")

        # Ensure nans are converted into -1
        image = np.nan_to_num(image, nan=-999.9)    
        # Create a mask based on the condition
        nan_mask3 = (image[:, :, 3] < 600)
        # Set invalid areas to -999.9 for all bands where nan_mask3 is True
        # image[nan_mask3] = -999.9
        # Initialize a validity mask with ones, same shape as the third channel of large_image
        validity_mask = np.ones_like(image[:, :, 3])
        # Set invalid areas (where nan_mask3 is True) to 0 in the validity mask
        validity_mask[nan_mask3] = 0
        # Expand the validity_mask to 3D by adding a new axis for the third dimension
        validity_mask = np.expand_dims(validity_mask, axis=-1)  # Shape will be (height, width, 1)
        #2019
        global_stats = [
            (1, 17632, 560.06472951343, 1196.0411021854),  # "Band_1"
            (1, 17008, 681.04916472418, 1198.6826622834),  # "Band_2"
            (1, 16448, 539.61406734262, 1205.572326046),   # "Band_3"
            (1, 15780, 1907.2003828845, 1511.6419398198),  # "Band_4"
            (1, 15226.5, 899.7007550269, 765.32225378592), # "Band_5"
            (1, 15097.5, 519.27680863559, 551.86427431662),# "Band_6"
            (1, 15910, 844.29700464677, 1205.2877346543),  # "Band_7"
            (1, 15881, 1621.8905432716, 1376.7183452063),  # "Band_8"
            (1, 15768.5, 1842.4491744618, 1463.592587156), # "Band_9"
            (1, 15627, 1968.4632191742, 1515.7663477273),  # "Band_10"
            (-49.931972503662, 21.146499633789, -22.631712949381, 10.497939369981), # "Band_11"
            (-39.526508331299, 30.792608261108, -13.370378733638, 7.5605403314878), # "Band_12"
            (-55.94380569458, 16.95357131958, -22.875294035681, 10.653140791668),   # "Band_13"
            (-42.266712188721, 29.136224746704, -14.138206542059, 8.4029188422164), # "Band_14"
            (126.5, 65535, 5656.8157776171, 3879.6824819856),  # "Band_15"
            (74.5, 57146, 2904.722165331, 2152.5637714749),    # "Band_16"
            (-27, 3077, 476.49736979174, 595.03618170533),     # "Band_17"
            (-129.1674041748, -122.98337554932, -125.29244682235, 1.515121082171), # "Band_18"
            (48.236663818359, 50.983592987061, 49.841552684935, 0.75116936146855), # "Band_19"
        ]

        
        # Apply normalization to each band
        for i in range(image.shape[2]):
            image[ :, :, i] = normalize_band(
                image[ :, :, i], 
                global_stats[i][2], 
                global_stats[i][3]
            )
        # Now concatenate along the last dimension
        image = np.concatenate((image, validity_mask), axis=-1)
        image = np.concatenate((image[ :, :, :4], image[ :, :, 10:14],image[ :, :, -1:]), axis=-1)
        # image = image[:, :, -5:]
        # Ensure nans are converted into -1
        image = np.nan_to_num(image, nan=-999.9)
        
        ###################################################################################
        #Predict using smooth blending
        
        # Use the algorithm. The `pred_func` is passed and will process all the image 8-fold by tiling small patches with overlap, 
        # called once with all those image as a batch outer dimension.
        # Note that model.predict(...) accepts a 4D tensor of shape (batch, x, y, nb_channels), such as a Keras model.
        predictions_smooth = predict_img_with_smooth_windowing(
            image,
            window_size=patch_size,
            subdivisions=2,  # Minimal amount of overlap for windowing. Must be an even number.
            nb_classes=2,#n_classes,
            pred_func=(
                lambda img_batch_subdiv: model.predict((img_batch_subdiv))
            )
        )
        
        print("prediction shape is {}".format(predictions_smooth.shape))
        
        meta = original_meta.copy()
        
        meta.update({
            'dtype': 'float32',
            'nodata': -999.9,
            # 'count': 2 if predictions_smooth.ndim == 2 else predictions_smooth.shape[2],
            'count': 2,
            'compress': 'lzw'
        })
        
        print("Appending...")
        all_preds.append(predictions_smooth)

    # Average predictions across all models
    avg_prediction = np.mean(all_preds, axis=0)
    # std_prediction = np.std(all_preds, axis=0)
    avg_prediction[nan_mask3, :] = -999.9
    # Create the prediction file path by appending the extracted number
    predicted_geotiff_path = os.path.join(output_pred_dir, f'prediction_{number}.tif')
    # # Write the predicted numpy array as a GeoTIFF
    with rasterio.open(predicted_geotiff_path, 'w', **meta) as dst:
        if avg_prediction.ndim == 1:
            dst.write(avg_prediction, 1)
        else:
            for i in range(avg_prediction.shape[2]):
                dst.write(avg_prediction[:, :, i], i + 1)
 
    print(f"Predicted GeoTIFF saved at: {predicted_geotiff_path}")
    # Clear any previous models or sessions
    

