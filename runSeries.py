#!/usr/bin/python

import os
import numpy as np

minfrag = []
numfrag = []
wCoh = []

os.system("rm output.txt")
reps = 100
for i in range(0,reps):
    print i+1,"of",reps
    os.system("rm results/*.log")
    os.system("./bin/main.exe >> output.txt")

with open("output.txt") as f:
    for line in f:
        if (line[0:20]=="Number of Fragments:"):
            l = line.strip().split()
            numfrag.append(float(l[-1]))
        if (line[0:21]=="Final Cohesive Energy"):
            l = line.strip().split()
            wCoh.append(float(l[-1]))
        if (line[0:22]=="Minimum fragment size:"):
            l = line.strip().split()
            minfrag.append(float(l[-1]))

assert(len(wCoh) == len(minfrag))
assert(len(wCoh) == len(numfrag))
#os.system("rm output.txt")
print ""
print "results of",len(minfrag),"experiments"
print "wCoh:    mean=",np.average(wCoh)," median=",np.median(wCoh)
print "numfrag: mean=",np.average(numfrag)," median=",np.median(numfrag)
print "minfrag: mean=",np.average(minfrag)," median=",np.median(minfrag)
