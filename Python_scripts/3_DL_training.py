#source /opt/conda/etc/profile.d/conda.sh conda -V
import os
import numpy as np
import random
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
from Deep_learning_UNET import get_unet_128_II#, get_unet_128
from loss_Wloss import rmse, bias, r_squared, laplace_loss#mae,


root_directory = '/home/barros/Big_tree/'
os.chdir(root_directory)

# Define function to extract unique identifier from the file name
import re
def extract_identifier(filename):
    match = re.search(r'\d+_\d+_\d+', filename)
    if match:
        return match.group(0)
    return None

# Visualize random image and mask pair from patches
train_images_dir = os.path.join(root_directory, '128_patches_NN/train/images')
train_masks_dir = os.path.join(root_directory, '128_patches_NN/train/masks')
val_images_dir = os.path.join(root_directory, '128_patches_NN/val/images')
val_masks_dir = os.path.join(root_directory, '128_patches_NN/val/masks')


# Get list of image and mask files
img_list = os.listdir(train_images_dir)
msk_list = os.listdir(train_masks_dir)

# Filter out files without a valid identifier
img_list = [f for f in img_list if extract_identifier(f) is not None]
msk_list = [f for f in msk_list if extract_identifier(f) is not None]

# Sort image and mask lists based on the extracted identifier from filenames
img_list = sorted(img_list, key=lambda x: extract_identifier(x))
msk_list = sorted(msk_list, key=lambda x: extract_identifier(x))

# Ensure the lists have the same number of elements
num_images = len(img_list)
num_masks = len(msk_list)
if num_images != num_masks:
    raise ValueError(f"Number of images ({num_images}) and masks ({num_masks}) do not match!")


# Visualize a random image and mask pair
img_num = random.randint(0, num_images - 1)

img_for_plot = tiff.imread(os.path.join(train_images_dir, img_list[img_num]))
mask_for_plot = tiff.imread(os.path.join(train_masks_dir, msk_list[img_num]))

image_name = extract_identifier(os.path.basename(img_list[img_num]))
mask_name = extract_identifier(os.path.basename(msk_list[img_num]))


image = tiff.imread(os.path.join(train_images_dir, img_list[img_num]))
image.shape

IMG_HEIGHT = image.shape[1]
IMG_WIDTH  = image.shape[2]
IMG_CHANNELS = image.shape[3]
n_classes=1

print("height= {}, Width = {}, Channels= {}, and classes = {}".format(IMG_HEIGHT,IMG_WIDTH,IMG_CHANNELS,n_classes))

print(img_for_plot.shape)
print(f"Selected image: {img_list[img_num]}")
print(f"Selected mask: {msk_list[img_num]}")
print(f"Image identifier: {image_name}")
print(f"Mask identifier: {mask_name}")


################################################################
# Define Generator for images and masks so we can read them directly from the drive. 

seed=24
num_class = 1
learning_rate=0.0001
decay_rate=0.1
batch_size = 32
activation_function = 'linear' #'sigmoid' is for classification with one class
num_train_imgs = len(os.listdir(train_images_dir))
num_val_images = len(os.listdir(val_images_dir))
epochs = 400
steps_per_epoch = num_train_imgs//batch_size
val_steps_per_epoch = num_val_images//batch_size

print(num_train_imgs,num_val_images,steps_per_epoch,val_steps_per_epoch)


def random_contrast(img, lower=0.95, upper=1.05):
    """Apply random contrast to the first 10 bands of the image."""
    # Check if the image has at least 10 bands
    if img.shape[3] < 10:
        raise ValueError("The input image must have at least 10 bands.")

    # Apply contrast adjustment only to the first 10 bands
    contrast_factor = np.random.uniform(lower, upper)
    mean = np.mean(img[:, :, :, :10], axis=(0, 1), keepdims=True)  # Calculate the mean of the first 10 bands

    # Apply contrast adjustment to the first 10 bands
    img[:, :, :, :10] = (img[:, :, :, :10] - mean) * contrast_factor + mean
    
    # Clip the values to the desired range
    img = np.clip(img, -3, 3)  # Adjust the clipping values based on your needs

    return img


# Custom data generator for segmentation tasks with data augmentation
def normalize_band(band, mean_val, std_val):
    # Normalize the band so that mean is 0, max is 1, and min is -1
    norm_band = (band - mean_val) / std_val  # Standardize
    return norm_band #np.clip(norm_band, -3, 3)  # Clip to ensure the values stay between -1 and 1

def min_max_scale_band(band, min_val, max_val):
     # Apply Min-Max scaling to the band
     scaled_band = (band - min_val) / (max_val - min_val)
     return np.clip(scaled_band, 0, 1)  # Ensure the values stay between 0 and 1

def preprocess_data(img, mask):
    img = np.nan_to_num(img, nan=-999.9)
    mask = np.nan_to_num(mask, nan=-999.9)

    # Create mask layer
    nan_mask1 = (mask <= 0.2)
    nan_mask2 = (mask > 100)
    nan_mask3 = (img[:, :, :, 3] < 450)

    mask[nan_mask1] = -999.9
    mask[nan_mask2] = 100
    mask[nan_mask3] = -999.9

    # Set invalid areas to -999.9 for all bands where nan_mask3 is True
    img[nan_mask3] = -999.9

    validity_mask = np.ones_like(mask)
    validity_mask[nan_mask1] = 0
    validity_mask[nan_mask2] = 0

    # # Slice the first 14 bands (before concatenating the mask)
    # img = img[:,:,:,:14]

    # Normalize each band based on the provided global statistics
    global_stats = [
        # Min, Max, Mean, Std, 99th Percentile for each band
        (1.0, 17632.0, 506.24448426415563, 1027.0768161896538, 10052.0),  # Band 1
        (1.0, 17008.0, 564.0686877027404, 1043.0635466938252, 10103.5),   # Band 2
        (1.0, 16448.0, 444.0549587648218, 1043.4553729753661, 9824.0),    # Band 3
        (1.0, 15780.0, 1407.1948212383818, 1503.5120523116918, 8248.0),   # Band 4
        (1.0, 15226.5, 687.836357071675, 728.7184260727269, 4355.5),      # Band 5
        (1.0, 15097.5, 414.29148988909753, 497.7401146494681, 3985.0),    # Band 6
        (1.0, 15910.0, 662.2172250632456, 1071.4429042615889, 9879.0),    # Band 7
        (1.0, 15881.0, 1207.5711119396537, 1340.5543037534756, 9217.5),   # Band 8
        (1.0, 15768.5, 1364.11075655541, 1452.1278275631232, 8514.0),     # Band 9
        (1.0, 15627.0, 1450.4100125633602, 1520.452667078374, 7818.0),    # Band 10
        (-76.40070343017578, 21.146499633789062, -26.954377849826436, 10.980594995524294, -3.003811244964614),  # Band 11
        (-46.556461334228516, 30.7926082611084, -15.11396322692858, 7.002571111634105, 7.6061601448059095),    # Band 12
        (-60.72349548339844, 16.953571319580078, -27.55455920900712, 11.124850667377768, -2.669262056350704),  # Band 13
        (-42.2667121887207, 29.1362247467041, -16.33617530904435, 7.746494639420921, 8.33240333557129)         # Band 14
    ]


    # Apply normalization to each band
    for i in range(img.shape[3]):
        img[:, :, :, i] = normalize_band(
            img[:, :, :, i], 
            global_stats[i][2], 
            global_stats[i][3]
        )

    # Concatenate validity mask as the last band
    img = np.concatenate((img, validity_mask), axis=-1)

    return img, mask


# Custom data generator for segmentation tasks with data augmentation
def custom_data_generator(img_dir, mask_dir, batch_size):
    img_files = sorted(os.listdir(img_dir))
    mask_files = sorted(os.listdir(mask_dir))

    # Data augmentation parameters
    img_data_gen_args = dict(horizontal_flip=True, vertical_flip=True,fill_mode='reflect')
    mask_data_gen_args = dict(horizontal_flip=True, vertical_flip=True,fill_mode='reflect')

    # Create ImageDataGenerators
    image_datagen = ImageDataGenerator(**img_data_gen_args)
    mask_datagen = ImageDataGenerator(**mask_data_gen_args)

    while True:
        for start in range(0, len(img_files), batch_size):
            end = min(start + batch_size, len(img_files))
            batch_img_files = img_files[start:end]
            batch_mask_files = mask_files[start:end]
            
            images = []
            masks = []

            for img_file, mask_file in zip(batch_img_files, batch_mask_files):
                img_path = os.path.join(img_dir, img_file)
                mask_path = os.path.join(mask_dir, mask_file)
                
                img = tiff.imread(img_path)
                mask = tiff.imread(mask_path)
                
                img, mask = preprocess_data(img, mask)

                    # Apply random brightness and contrast augmentations
                rn = np.random.uniform()
                if rn > 0.75:
                    img = random_contrast(img)           
                # Ensure the dimensions are (128, 128, 16) and (128, 128, 1)
                if img.ndim == 4 and img.shape[0] == 1:
                    img = img[0]
                if mask.ndim == 4 and mask.shape[0] == 1:
                    mask = mask[0]
                angle = np.random.choice([0,-90, 90])
                if angle != 0:
                    img = rotate(img, angle, axes=(0, 1), reshape=False, mode='reflect')
                    mask = rotate(mask, angle, axes=(0, 1), reshape=False, mode='reflect')

                # Ensure last img band is binary and mask has no negative values
                img[:,:,-1] = np.where(img[:,:,-1] >= 0.5, 1.0, 0.0).astype(np.float32)
                mask = np.where(mask >= 0, mask, -999.9).astype(np.float32)

                # Ensure that there are no nan valeus
                img = np.nan_to_num(img, nan=-999.9)
                mask = np.nan_to_num(mask, nan=-999.9)
                
                images.append(img)
                masks.append(mask)
            
            images = np.array(images)
            masks = np.array(masks)
            
            # Apply data augmentation
            seed = np.random.randint(1e6)
            aug_images = image_datagen.flow(images, batch_size=batch_size, seed=seed, shuffle=False)
            aug_masks = mask_datagen.flow(masks, batch_size=batch_size, seed=seed, shuffle=False)
            
            yield next(aug_images), next(aug_masks)


# Set path to training and validation data
train_img_path = os.path.join(root_directory, '128_patches_NN/train/images')
train_mask_path = os.path.join(root_directory, '128_patches_NN/train/masks')


val_img_path = os.path.join(root_directory, '128_patches_NN/val/images')
val_mask_path = os.path.join(root_directory, '128_patches_NN/val/masks')


# Hyperparameters
batch_size = 32
seed = 24
replace_value = -999.9

# Create data generators
train_img_gen = custom_data_generator(train_img_path, train_mask_path, batch_size)
val_img_gen = custom_data_generator(val_img_path, val_mask_path, batch_size)
test_img_gen = custom_data_generator(test_img_path, test_mask_path, batch_size)




###########################################################################
#Define the model metrics and load model. 
metrics=[laplace_loss,rmse, bias, r_squared] 


# Initial and minimum learning rates
learning_rate = 0.0001  # ηmax
min_learning_rate = 0  # ηmin

# Calculate total training steps
total_steps = epochs * steps_per_epoch*0.95

# Define the cosine annealing learning rate schedule
lr_schedule = CosineDecay(
    initial_learning_rate=learning_rate,  # ηmax
    decay_steps=total_steps,  # Tmax
    alpha=min_learning_rate / learning_rate  # α = ηmin / ηmax
)


optimizer = Adam(learning_rate=lr_schedule)


IMG_HEIGHT = image.shape[1]
IMG_WIDTH  = image.shape[2]
IMG_CHANNELS = image.shape[3]+1
n_classes=1

print("height= {}, Width = {}, Channels= {}, and classes = {}".format(IMG_HEIGHT,IMG_WIDTH,IMG_CHANNELS,n_classes))




# Compile the model
def get_model():
    return get_unet_128(n_classes=n_classes, IMG_HEIGHT=IMG_HEIGHT, IMG_WIDTH=IMG_WIDTH, IMG_CHANNELS=IMG_CHANNELS, activation = activation_function)

model = get_model()

model.compile(optimizer=optimizer,
                  loss=laplace_loss, 
              metrics= metrics)

model.summary()



# Set up callbacks
callbacks_list = [
    CSVLogger("Model_SEN1_2/epoch_history.csv", separator=";", append=False),
    ModelCheckpoint(
        filepath="Model_SEN1_2/unet_tf2_{epoch:03d}_{val_laplace_loss:.4f}.h5",
        monitor='val_laplace_loss',  # Monitor the validation F1 score
        save_best_only=True,
        save_weights_only=True,
        mode="min",  # F1 score needs to be maximized
    ),
    EarlyStopping(
        monitor='val_laplace_loss',
        patience=400,
        mode="min"  # F1 score needs to be maximized
    )
]



# # Model training with augmentation

# Fit the model
history1 = model.fit(
    train_img_gen,
    epochs=epochs,
    steps_per_epoch=steps_per_epoch,
    validation_data=val_img_gen,
    validation_steps=val_steps_per_epoch,  # Corrected parameter name
    callbacks=callbacks_list,  # Callbacks should be a list of callback functions
    shuffle=False,  # Optionally set to True if you want to shuffle the data
    verbose=1  # Set verbosity level (0, 1, or 2)
)



print('This is the end of the script')
