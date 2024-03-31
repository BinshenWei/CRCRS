# Exploring prognostic biomarkers in pathological images of colorectal cancer patients via deep learning
### Abstract
- Background：Whole-slide images (WSIs) of cancer patients provide valuable information for disease progression and patient survival. However, extracting prognostic indicators from pathology images is challenging due to the subtle and complex nature of phenotypic information.
- Methods：In this study, we developed a weakly supervised deep learning model using a cohort of 640 colorectal cancer (CRC) patients from the PLCO dataset and 522 CRC patients from the TCGA dataset. Our goal was to identify new prognostic markers for pathological WSIs. We created the Colorectal Cancer Risk Score (CRCRS) to assess patient prognosis and visualized the pathological phenotype of CRCRS using Grad-CAM. Additionally, we utilized multiomics data from the TCGA CRC cohort to investigate the potential pathogenesis of CRCRS.
- Results：Our survival analysis revealed that CRCRS was an independent prognostic factor for both the PLCO (p<0.001) and TCGA (p<0.001) cohorts, and its predictive power was independent of the clinical staging system. Combining CRCRS with the TNM staging systems improved the accuracy of the nomogram in predicting patient prognosis compared to using the TNM staging system alone. We identified features of CRCRS, including immature tumor mesenchyme, disordered myxoid gland structures, clusters of small tumor cells, and infiltrating inflammatory cells, as the main underlying cells. Analysis of TCGA CRC data suggested associations between CRCRS and the activation of energy production and metabolic pathways, the tumor immune microenvironment, and genetic alterations in APC, SMAD2, EF1AKMT4, and EPG5.
- Conclusion：Our deep learning model established the CRCRS as a prognostic predictor for CRC patients. It provides a valuable method for CRC risk stratification and precision treatment of patients.

![Fig 1](https://github.com/BinshenWei/CRCRS/assets/162391399/98a3af63-37d1-42ec-bfb4-1a8dcf024999)
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
TiatoolBox 1.5.1
DigiPath 0.1.5

```
### Data Preparation:

We can download the H&E slides from the TCGA and PLCO datasets. 

The TCGA H&E WSIs can be downloaded in https://portal.gdc.cancer.gov/. You can use the cdc-client tool to download these H&E slides. 

The PLCO H&E WSIs can be downloaded in https://cdas.cancer.gov/plco/. You need to have request to access  this dataset. 

The TCGA survival information can be downloaded in https://www.cbioportal.org/ and the PLCO survival information can be downloaded in https://cdas.cancer.gov/plco/. 



### Tumor region segmentation network

In this study, we employed a pretrained U-Net tissue segmentation model (https://github.com/haranrk/DigiPathAI/blob/master/DigiPathAI/Segmentation.py) to delineate tumor regions (Fig. 1a )The U-Net tumor segmentation model was trained on the DigestPath-2019 dataset, which includes 872 images of CRC tissue. Our research involved the use of the U-Net algorithm to identify tumor regions in WSIs via a sliding window approach. Specifically, we determined the sliding window dimensions to be 1024 x 1024 (40x magnification) and 512 x 512 (20x magnification). Subsequently, we adjusted the size of the patches from 1024 x 1024 to 512 x 512 and imported them into the U-Net model.

For each slide, we calculate the tumor ratio (tumor area / patch area) of each patch. Then, we sort these patches by using the tumor ratios and select the top k patches to train the CRC prognostic network. 

### Dataloader 

The detail can be seen in <code>survival.py</code> . We select the top k tumor ratio patch. And we used the ResNet backbone to extract features of each patch.  Finally, we get the features of each slide, overall survival information. 

### Prognostic network

We use the TransMIL to extract features and used the Cox regression to calculate the loss. In this study, we assign the k to {4, 8, 16, 32} and test the predictive performance with different k. 

We can see the detail in the model from the <code>Prognostic.py</code>



### Prognostic network Training 

You can set the different parameters by your actual situations, such as the segmentation and H&E slides path. Then, you can use <code>python train.py</code>. 

### Prognostic network valid

You can set the different parameters by your actual situations, such as the model checkpoint path. Then, you can use <code>python test.py</code>. 

### Script outputs:
.xlsx file with CRCRS(Colorectal Cancer Risk Score), patients' OS, outcome, patient ID.
To generate heatmap , use script . for help regarding input parameters.



### Statistical analysis

See the Statistical analysis folder. 
