from sklearn.ensemble import StackingClassifier
from sklearn.svm import SVC
from sklearn.ensemble import (RandomForestClassifier, AdaBoostClassifier,GradientBoostingClassifier, ExtraTreesClassifier)
estimators = [('rf_model',RandomForestClassifier()),
              ('ab_model', AdaBoostClassifier()),
              ('ex_model',ExtraTreesClassifier()),
              ('gdbc_model',GradientBoostingClassifier()), 
              ('svc_model',SVC())
]
clf = StackingClassifier(
    estimators=estimators, final_estimator=SVC(probability=True),cv=10
)
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
    new_data, Y, test_size=0.2)
clf.fit(X_train, y_train)
score=clf.score(X_test, y_test)
print(score)
