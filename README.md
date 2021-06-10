# PPDTS
Predicting drug-target interactions based on network similarity

# Requirements
* python 3.6.12
* numpy 1.19.2
* pandas 1.1.1 
* scipy 1.5.2
* sklearn 
* pytorch 1.4
* matplotlib 3.3.1

# Data description--two datasets

* data set 1 (Drugs: 786; Targets: 527; Associations:4547). 
* data set 2 (Drugs: 575; Targets: 981; Associations:8008).
* MATADOR_chem_protein_direct.xlsx: 8936 drug-target interactions extracted from MATADOR.

# Code description

* The ppdts.py: Code representation of PPDTS function implementation process.

* The k_cv_5.py: 5-k cross validation method involved in the experiment.

* The calcul_height_score.py: The Target set to obtain the highest score of the Drug set in the predicted results by PPDTS.

# results description

* 786_predictd_result.xlsx: From the data set 1, PPDTS predicted results.
* new_37_DTIs.xlsx: Among 786 DTIs, 37 were verified by UniProt and DrugBank databases.

# Supplementary Information description
* Supplementary Information.docx: Figure S1: the proportion of negative samples in the total category; Figure S2: different combination weights α corresponding PPDTS experimental results; Figure S3: AUC curve corresponding to data set 1;
Table S1: experimental evaluation parameters corresponding to dataset 2.

# If you have any questions, please contact： 

* wangyongqing422@gmail.com
