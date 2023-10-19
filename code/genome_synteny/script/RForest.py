import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import log_loss
from random import shuffle

np.random.seed(0)

X = []
y = []
strline = []
infile = open('out.txt','rb')
for line in infile:
    strline.append(line)

infile.close()
shuffle(strline)
for line in strline:
    line = line.strip('\n')
    sent = line.split('\t')
    xi = sent[2:14]
    yi = int(sent[14])
    xi = np.asarray(xi,dtype = float)
    X.append(xi)
    y.append(yi)

X = np.asarray(X)
y = np.asarray(y)

X_train, y_train = X[:12000], y[:12000]
X_valid, y_valid = X[12000:16000], y[12000:16000]
X_train_valid, y_train_valid = X[:16000], y[:16000]
X_test, y_test = X[16000:], y[16000:]

clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train_valid, y_train_valid)
clf_probs = clf.predict_proba(X_test)
score = log_loss(y_test, clf_probs)
print clf.score(X_test, y_test)
rdata = []
infile = open('imdata.txt','rb')
for line in infile:
    line = line.strip('\n')
    sent = line.split('\t')
    xi = sent[2:14]
    xi = np.asarray(xi,dtype = float)
    rdata.append(xi)

rdata = np.asarray(rdata)
yy_pre = clf.predict(rdata)
yy_pro = clf.predict_proba(rdata)
np.savetxt("result.txt",yy_pre)
np.savetxt("result_pro.txt",yy_pro)
np.savetxt("important_feature.txt",clf.feature_importances_)
