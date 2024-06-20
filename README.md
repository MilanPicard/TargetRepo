# TargetRepo
Target repositioning on advanced prostate cancer using a multi-layer network and machine learning.
This is the code and data necessary to reproduce the results presented in the following article.

### To cite
Picard, Milan, et al. "Target repositioning using multi-layer networks and machine learning: the case of prostate cancer." Computational and Structural Biotechnology Journal (2024). https://doi.org/10.1016/j.csbj.2024.06.012


### Available functions
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
