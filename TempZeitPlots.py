# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:21:25 2022

@author: Bernhard
"""
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

s1=500
s2=450

b1=50
b2=50

k1=0.6
k2=0.4

def f1(x):
    return s1-(s1-b1)*math.exp(-k1*x)

def f2(x):
    return s2-(s2-b2)*math.exp(-k2*x)

xval=np.arange(0,8,0.05)

y1=np.array(list(map(f1, xval)))
y2=np.array(list(map(f2, xval)))

fig1 = plt.figure(figsize=(7,5))
ax1  = fig1.add_subplot(1,1,1)

ax1.plot(xval, y1, 'k', label="ungekühlt")
ax1.plot(xval, y2, 'b', label="gekühlt")
legend = ax1.legend(loc="lower right", fontsize="x-large")
ax1.yaxis.set_major_formatter(plt.NullFormatter())
plt.xlabel("Zeit")