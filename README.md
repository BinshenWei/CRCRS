[Fig1.pdf](https://github.com/BinshenWei/CRCRS/files/14588132/Fig1.pdf)
# Exploring prognostic biomarkers in pathological images of colorectal cancer patients via deep learning
### Abstract
- Background
Whole-slide images (WSIs) of cancer patients offer significant insights into disease advancement and patient outcome prediction. Nevertheless, the extraction of prognostic markers from pathology images is intricate due to the nuanced and intricate characteristics of phenotypic data.
- Methods
In this research, a weakly supervised deep learning model was formulated utilizing a group of 640 colorectal cancer (CRC) patients sourced from the PLCO dataset and 522 CRC patients from the TCGA dataset. The primary objective was to discover novel prognostic indicators for pathological WSIs in CRC. A novel metric termed the Colorectal Cancer Risk Score (CRCRS) was devised to evaluate patient prognosis, and the pathological phenotype of CRCRS was visualized using Grad-CAM. Furthermore, multi-omics data from the TCGA CRC cohort was leveraged to explore the potential pathogenesis associated with CRCRS.
- Results
The results of our survival analysis indicated that CRCRS served as an independent prognostic indicator for both the PLCO (p<0.001) and TCGA (p<0.001) cohorts, with its predictive efficacy remaining unaffected by the clinical staging system. Integration of CRCRS with the TNM staging systems resulted in enhanced accuracy in predicting patient prognosis compared to utilizing the TNM staging system in isolation. Noteworthy features of CRCRS, such as immature tumor mesenchyme, disorganized myxoid gland structures, small clusters of cancer cells, and infiltrating inflammatory cells, were identified. Examination of TCGA CRC multi-omics data revealed potential correlations between CRCRS and the activation of energy production and metabolic pathways, the tumor immune microenvironment, as well as genetic mutations in APC, SMAD2, EF1AKMT4, EPG5 and TANC1.
- Conclusion
Our deep learning algorithm identified CRCRS as a prognostic indicator for patients with colorectal cancer. This finding offers a significant approach for stratifying the risk of colorectal cancer and tailoring precise treatment strategies for individual patients.



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
TIAToolbox 0.1.5
Digipath 0.1.5

```
### Classification network:
#### Data Preprocess:
Run python ./data_preprocess/cls/make_slide_cuting.py to generate the thumbnail.
Run python ./data_preprocess/cls/make_tissue_mask.py to generate the tissue mask.
Run python ./data_preprocess/cls/make_lab_mask.py and python ./data_preprocess/cls/make_other_mask.py to generate Tumor,
Liver, Stroma and Hemorrhage & Necrosis mask.
Run python ./data_preprocess/cls/create_patch.py to obtain patches.

Training data should be a dictionary stored to disk containing the following keys set in ./configs/resnet18_tcga.json:
"data_path_40": the path to store tiled slides in non-overlapping 768x768 pixel windows at a magnification of 40x
"json_path_train":the path to store the training json containing the dictionary of a list of tuple (x,y) coordinates. 
"json_path_valid":the path to store the valid json containing the dictionary of a list of tuple (x,y) coordinates. 
Valid data should contain the slides and masks of the foreground area for each slide.

#### Classification network :
To train a model, use script ./Classification/bin/train_tcga.py. Run python ./Classification/bin/train_tcga.py -h to get help regarding input parameters.
The data path and label path are set in ./configs/resnet18_tcga.json (data_path_40,  json_path_train, json_path_valid).

#### Prognostic network training:
To train the prognostic network, use script ./Prognositic/bin/train.py. Run python  ./Prognositic/bin/train.py -h for help regarding input parameters.

#### Prognostic network valid:
To run a model on a test set, use script ./Prognositic/vis/vis_save.py. Run python  ./Prognositic/vis/vis_save.py -h for help regarding

#### Script outputs:
.csv file with TRS(tumor risk score), patients' OS, outcome, patient ID.
To generate heatmap , use script ./Prognositic/vis/vis_save_1.py. Run python  ./Prognositic/vis/vis_save_1.py -h for help regarding input parameters.



