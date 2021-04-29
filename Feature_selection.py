import numpy as np
from sklearn.linear_model import Lasso,LassoCV
model_lasso=LassoCV(cv=5, alphas=np.array([0.001,0.002,0.003,0.004, 0.005, 0.006, 0.007, 0.008,0.009, 0.01])).fit(X, Y)
model_lasso.fit(X,Y)
mask = model_lasso.coef_ != 0 
new_data = X[:,mask]
