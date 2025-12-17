# The Last Giants: Deep Learning Reveals Critical Conservation Gaps in Canada's Coastal Temperate Rainforests 

**Authors:**\
Luizmar de Assis Barros<sup>a</sup>, Karen Price<sup>b</sup>, José Bermúdez<sup>c</sup>, Chris Johnson<sup>a</sup>, Camile Sothe<sup>d</sup>, Juan Pablo Ramírez-Delgado<sup>a</sup>, Xavier Llano<sup>a</sup>, Alemu Gonsamo<sup>b</sup>, Michelle Venter<sup>a</sup>, Oscar Venter<sup>a</sup>

a *University of Northern British Columbia, 3333 University Way, Prince George, V2N 4Z9, British Columbia, Canada*\
b *Independent Researcher, Salt Spring Island, V8K 1Y8, Canada*\
c *School of Earth, Environment & Society, McMaster University, Hamilton, L8S 4K1, Canada*\
d *Planet Labs PBC, San Francisco, 695571, California, USA*

---

## Abstract

Canada’s coastal temperate rainforests support some of the world’s largest trees. Yet, decades of logging concentrated in the most productive forests, combined with the protection of less-productive areas, have transformed much of this ecosystem into younger, homogeneous stands. Such biases, compounded by natural mortality and slow recruitment, place big trees among the most imperiled organisms on Earth. Effective conservation efforts must address these historical biases by protecting rare big-treed forests. Current inventories, while valuable, lack spatial resolution, accuracy, and the ability to capture forests sustaining the largest trees. To overcome these limitations, we developed a high-resolution (10-m) big-treed forest map by integrating airborne LiDAR, satellite imagery, and a deep learning approach. From the 2.7 million hectares of forests, we identified approximately 133,000 hectares of big-treed forests. Roughly half occurred within human-modified landscapes and half within regions with reduced accessibility and limited disturbance. Of these last remaining big-treed forests, 62% were unprotected. Notably, less than a quarter of the existing protected areas (PAs) and 27% Other Effective area-based Conservation Measures (OECMs) contain at least 1 ha of big-treed forest. Together, these results show that big-treed forests are under-protected, highly fragmented, and vulnerable to further human-caused loss. Our study provides a scalable, open-access framework for high-accuracy mapping and monitoring of big-treed forests to support targeted conservation of these rare and vulnerable ecosystems.

---

## Significance Statement

Big-treed forests are among the most ecologically valuable yet most threatened components of coastal temperate rainforests. Their conservation is hindered by the lack of accurate, high-resolution maps capable of identifying areas where the largest trees persist and thus inform policy. We address this gap by developing a 10-m big-treed forest map using airborne LiDAR, satellite imagery, and deep learning, revealing their true extent, distribution, and protection status. Our results show that these forests are rare, highly fragmented, and disproportionately unprotected, with about half occurring in remote regions with reduced-access and half in human-modified landscapes. Although our study focused on Vancouver Island, Canada, our approach provides a scalable, open-access framework to guide conservation planning and long-term monitoring of threatened big-treed ecosystems worldwide.

---

## Keywords
Keystones, LiDAR, Multispectral Imagery, Old-growth forests, SAR imagery.

---

## Main Result

![Big-tree Forests of Vancouver Island](Main_fig3.jpg)
**Figure 1** Bivariate map of big tree forests locations, classified into top 10% (> 44m), top 5% (>=48m) and 1% (>=54m) tallest canopies and old (>=250 years); mature (>80 and <250 years); and young/second growth forest (<80years) on the Islands on Vancouver Island: a) large unprotected old big-tree forests northwest of the island; b) old-growth management areas (OGMAs), a type of OECM, north of the Woss village in the Nimpkish Valley; c) southeast Strathcona park boundary; d) Coastal Douglas fir mature big-tree forest inside the Saysutshun (Newcastle Island Marine) Park; and e) Fairy Creek watershed partially covered by OGMAs. 

# Python Scripts Description:
**- 1_Prep_IMG_MSK.ipynb**\
  Download and pre-processing of the wall-to-wall data predictors, as well as alignment all predictors and reference data tiles\
**- 2_Patchify_IGM_MSK.ipynb**\
  Sample acquisition of coregistered patches (128x128 pixels, 10m resolution) of the reference and predictors data\
**- 3_DL_training.py**\
  Deep learning model training using training and validation patches obtained from step 2\
**- 4_Deep_learning_UNET.py**\
  Two U-Net architectures used during our ablation study. While both structures are nearly identical, the second includes a second output channel for uncertainty estimation\
**- 5_loss_Wloss.py**\
Loss functions used for model training and assessment\
**- 6_utils.py**\
Functions used for reading and processing image patches, obtained from https://github.com/Vooban\
**- 7_Predict_height.py**\
Once the model was calibrated and the best five epochs selected, we predicted canopy height and cover using our wall-to-wall predictors. We loaded each of the best five epochs, predicted, and combined the results.\
**- 8_smooth_tile_predition.py**\
Predictions results often degraded towards the edges of patches. This function helps minimize this effect for a seamless merging of all predicted patches. Obtained fro https://github.com/Vooban\
**- 9_Crop&merge_predition.ipynb**\
To further improve prediction on edges, our prediction was performed for tiles with 250-pixel buffers. Here, we removed the buffer and merged the results.

# R Scripts Description:
**- 1_Big_treed_Forest_Maps.R**\
  Here we mapped big-treed forest using our locally trained canopy cover  and height layers, as well as an alternative global dataset and a local map of forest inventory commonly used in local studies\
**- 2_Fig1_FigS1_S2_S6_S8.R**\
  Data processing leading to Figure 1 and Supplemental Figures 1, 2, 6 and 8. Here we obtain ~1.1 million reference LiDAr pixels and use to assess models' performance and residuals\
**- 3_FigS3_Tables.R**\
  Processing of all big-treed forest maps leading to Figure S3 and main and supplemental result tables. Here we assessed big-treed forest overlaps across different scales and alternative dataset\
**- 4_FigureS4_S5.R**\
  Processing of big-treed forest centroids to assess big-treed forest distribution across environmental gradients. 

