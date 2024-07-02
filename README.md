# Early-prediction-of-Antimicrobial-Resistance-AMR-using-Machine-Learning

# Summary
Using genetic information from resistant bacteria, it is possible to model the characteristics that confer resistance to certain antibiotics. This information could be used to predict new potential resistant genes or sets of genes and take appropriate action. In this project, real data from a series of resistant bacteria identified in the laboratory will be used to build a machine learning model capable of predicting potential resistance genes (ARGS).

# Repository structure
data/: contains input data for the project (not all data is included in the repository due to size limitations)
eda/: Exploratory Data Analysis (EDA) of the input data
data_cleaning: data cleaning and preparation of the input data for the ML models
correlation_analysis: correlation analysis between the different variables in the dataset
ml/: contains all the Machine Learning (ML) related work
ML_for_AMR_prediction: machine learning training and evaluation of different models to predict Antimicrobial Resistance (AMR)
workflows/: contains ETL workflows to obtain the input data for the project
card: identifies ARGs and SNPs present in the samples by CARD.
data_collection_bvbrc/: contains various scripts to retrieve the input data from BV-BRC source
data_collection_ncbi/: contains various scripts to retrieve the input data from NCBI source
resfinder/: identifies the ARGs present in the samples by Resfinder v4.3.1
config.yaml: configuration file for the workflows
