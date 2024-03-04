# TargetRepo
Target repositioning on advanced prostate cancer using a multi-layer network and machine learning

utils.r
general functions to run the code

extract_functions.R
functions used to extract features from the graph

feature_extraction_functions.Rmd
Code that extracted features from the graph

Feature_selection.Rmd
Split the data into training and test dataset, then feature selection on the training set

Machine_learning.Rmd
Training and tuning of 5 machin learning models, evaluation on the test set

Final_prediction.Rmd
Using the best model (KNN) and the best features (Non directional random walks and topological similarities), 
train the model on all data and predict new potential targets
