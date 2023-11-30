# ALS Survival Prediction
This is a repository of the code used for the paper ["Accurate Personalized Survival Prediction for Amyotrophic Lateral Sclerosis Patients"](https://www.nature.com/articles/s41598-023-47935-7) by Kuan et al. We will upload the codebase once the paper is published. Note that the data will not be uploaded.

This repository is modified from [ISDEvaluation](https://github.com/haiderstats/ISDEvaluation) by Humza Haider.

# How to make ISD prediction from new patient using our pre-triained model?
1. Put the patients' input variables in the "pretrained/new_input.csv" file—one row for a patient. Refer to the "Data_lagend.xlsx" for the variable definition. (Note that the "new_input.csv" file only contains the variables required by our final pre-trained model. We consider all variables, including MR image, when building our the pre-trained model. However, the software decides not to include some variables because they do not improve the model performance based on its statistical results. Please refer to the paper for more description.)

2. Run the "new_patients.R" script. The script will plot the survival curves for every patient (in one plot) and write the "pretrained/result_curves.csv" file. The first column of the "result_curves.csv" is the time point, and after that, one column for a patient; each row is the survival probability of the corresponding time point.