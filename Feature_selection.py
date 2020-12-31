import numpy as np
from sklearn.linear_model import ElasticNet,ElasticNetCV
enet = ElasticNetCV(alphas=np.array([0.01, 0.02, 0.03,0.04, 0.05, 0.06, 0.07, 0.08,0.09, 0.1]),l1_ratio=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]).fit(X,Y)
mask = enet.coef_ != 0
new_data = X[:,mask]
print(new_data.shape)
