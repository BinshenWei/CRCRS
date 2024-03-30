# Exploring prognostic biomarkers in pathological images of colorectal cancer patients via deep learning
### Abstract
- Background：Whole-slide images (WSIs) of cancer patients provide valuable information for disease progression and patient survival. However, extracting prognostic indicators from pathology images is challenging due to the subtle and complex nature of phenotypic information.
- Methods：In this study, we developed a weakly supervised deep learning model using a cohort of 640 colorectal cancer (CRC) patients from the PLCO dataset and 522 CRC patients from the TCGA dataset. Our goal was to identify new prognostic markers for pathological WSIs. We created the Colorectal Cancer Risk Score (CRCRS) to assess patient prognosis and visualized the pathological phenotype of CRCRS using Grad-CAM. Additionally, we utilized multiomics data from the TCGA CRC cohort to investigate the potential pathogenesis of CRCRS.
- Results：Our survival analysis revealed that CRCRS was an independent prognostic factor for both the PLCO (p<0.001) and TCGA (p<0.001) cohorts, and its predictive power was independent of the clinical staging system. Combining CRCRS with the TNM staging systems improved the accuracy of the nomogram in predicting patient prognosis compared to using the TNM staging system alone. We identified features of CRCRS, including immature tumor mesenchyme, disordered myxoid gland structures, clusters of small tumor cells, and infiltrating inflammatory cells, as the main underlying cells. Analysis of TCGA CRC data suggested associations between CRCRS and the activation of energy production and metabolic pathways, the tumor immune microenvironment, and genetic alterations in APC, SMAD2, EF1AKMT4, and EPG5.
- Conclusion：Our deep learning model established the CRCRS as a prognostic predictor for CRC patients. It provides a valuable method for CRC risk stratification and precision treatment of patients.
### Prerequisites:
```
Python 3.7
Numpy 1.17.2
Scipy 1.3.0 
Pytorch 1.3.0/CUDA 10.1
torchvision 0.4.1
Pillow 6.0.0
opencv-python 4.1.0.25
openslide-python 1.1.1
Tensorflow 2.0.0
Tensorboard 2.0.0
tensorboardX 1.7
pandas 0.25.3
lifelines 0.23.4
### Data Preparation:
Use gdc_manifest.2020-01-03.txt to download WSI from TCGA dataset.
Use the download tool 'gdc-client' to download data. Run the command:
gdc-client download -m gdc_manifest.2020-01-03.txt -d ./data/svs/
### Classification network:
#### Data Preprocess:
Run python ./data_preprocess/cls/make_slide_cuting.py to generate the thumbnail.
Run python ./data_preprocess/cls/make_tissue_mask.py to generate the tissue mask.
Run python ./data_preprocess/cls/make_lab_mask.py and python ./data_preprocess/cls/make_other_mask.py to generate Tumor,
Liver, Stroma and Hemorrhage & Necrosis mask.
Run python ./data_preprocess/cls/create_patch.py to obtain patches.




