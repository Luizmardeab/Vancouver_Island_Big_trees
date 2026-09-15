# Local and Global Canopy Height Models Reveal Consistent Large-Tree Forest Conservation Gaps Despite Low Spatial Agreement

**Authors:**\
Luizmar de Assis Barros<sup>a</sup>, Karen Price<sup>b</sup>, Chris Johnson<sup>a</sup>, José Bermúdez<sup>c</sup>, Juan Pablo Ramírez-Delgado<sup>a</sup>, Xavier Llano<sup>a</sup>, Camile Sothe<sup>d</sup>, Alemu Gonsamo<sup>b</sup>, Michelle Venter<sup>a</sup>, Oscar Venter<sup>a</sup>

a *University of Northern British Columbia, 3333 University Way, Prince George, V2N 4Z9, British Columbia, Canada*\
b *Independent Researcher, Salt Spring Island, V8K 1Y8, Canada*\
c *School of Earth, Environment & Society, McMaster University, Hamilton, L8S 4K1, Canada*\
d *Planet Labs PBC, San Francisco, 695571, California, USA*

---

## Abstract

Canada’s coastal temperate rainforests are home to some of the world’s largest trees. These forests support disproportionate carbon storage, high structural complexity, and habitat for a wide range of species. Historical logging focused on the most productive stands, combined with the protection of less productive areas, has transformed much of this ecosystem into younger, more homogeneous forest landscapes. Effective conservation therefore requires identifying and protecting remaining large-tree forests. Medium- and high-resolution canopy height products offer a practical proxy for detecting these ecosystems. Here, we evaluated the ability of locally calibrated and global canopy height products to detect forests supporting the largest trees. We proposed an approach for identifying large-tree forests from canopy height models and applied it across ~2.7 million ha of forests on Vancouver Island, British Columbia, Canada. We identified approximately 133,000 ha of large-tree forests, with nearly half occurring within highly modified landscapes. Spatial agreement among local and global canopy height models was low (16–47%), indicating substantial disagreement in the mapped location of large-tree forests. Despite these differences, all models consistently revealed major conservation gaps, with 69–76% of large-tree forests occurring outside formal protected areas. The protected portion of these forests was highly fragmented, with fewer than 25% of protected areas containing more than 1 ha of large-tree forest. Together, these results show that large-tree forests are under-protected, fragmented, and vulnerable to human-caused loss. Our study provides an operational framework for identifying and monitoring these rare and vulnerable ecosystems at broad spatial scales.

---

## Keywords
Airborne LiDAR, Big Trees, Deep Learning, Forest Conservation, Other Effective area -based conservation measures (OECMs), Old-growth Forests, Synthetic Aperture Radar.

---

## Main Result

![Big-tree Forests of Vancouver Island](Main_fig3.jpg)
**Figure 1** Bivariate map of big tree forests locations, classified into top 10% (> 44m), top 5% (>=48m) and 1% (>=54m) tallest canopies and old (>=250 years); mature (>80 and <250 years); and young/second growth forest (<80years) on the Islands on Vancouver Island: a) large unprotected old big-tree forests northwest of the island; b) old-growth management areas (OGMAs), a type of OECM, north of the Woss village in the Nimpkish Valley; c) southeast Strathcona park boundary; d) Coastal Douglas fir mature big-tree forest inside the Saysutshun (Newcastle Island Marine) Park; and e) Fairy Creek watershed partially covered by OGMAs. 

# Data availability:
Data generated in this study are archived on Zenodo (DOI:10.5281/zenodo.17992142 or <https://zenodo.org/records/17992142>)

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
Functions used for reading and processing image patches, obtained from <https://github.com/Vooban>\
**- 7_Predict_height.py**\
Once the model was calibrated and the best five epochs selected, we predicted canopy height and cover using our wall-to-wall predictors. We loaded each of the best five epochs, predicted, and combined the results.\
**- 8_smooth_tile_predition.py**\
Predictions results often degraded towards the edges of patches. This function helps minimize this effect for a seamless merging of all predicted patches. Obtained from <https://github.com/Vooban>\
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

