import numpy as np
import os, glob

logdir = "experiments/log"
X1, X2, X3 = [], [], []

for log in glob.glob(os.path.join(logdir, "1.[0-9]*.txt")):
    logtext = open(log).read().split('\n')[:-1]
    logtext = [x.split() for x in logtext]
    X1.append(np.matrix(logtext, dtype="float32"))

for log in glob.glob(os.path.join(logdir, "2.[0-9]*.txt")):
    logtext = open(log).read().split('\n')[:-1]
    logtext = [x.split() for x in logtext]
    X2.append(np.matrix(logtext, dtype="float32"))

for log in glob.glob(os.path.join(logdir, "3.[0-9]*.txt")):
    logtext = open(log).read().split('\n')[:-1]
    logtext = [x.split() for x in logtext]
    X3.append(np.matrix(logtext, dtype="float32"))


X1_avg = np.mean(X1, axis=0)[:, ::-1]
X2_avg = np.mean(X2, axis=0)[:, ::-1]
X3_avg = np.mean(X3, axis=0)[:, ::-1]


for x in X1_avg:
    print x

print
for x in X2_avg:
    print x

print
for x in X3_avg:
    print x




