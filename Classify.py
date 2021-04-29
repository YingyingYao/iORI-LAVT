from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import VotingClassifier
clf1=GaussianNB()
clf2=LogisticRegression()
vote_clf=VotingClassifier(estimators = [('NB',clf1),('LOG',clf2)],voting='soft')
for clf in (clf1,clf2,vote_clf):
    score=cross_val_score(clf,X,Y,cv=10).mean()
    print(clf.__class__.__name__,score)
