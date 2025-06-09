# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:21:47 2024

@author: Xavier Guidetti
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np

# Getting back the objects:
with open('dists0_5.pkl','rb') as f:  # Python 3: open(..., 'rb')
    dist05 = pickle.load(f)
    
with open('dists5.pkl','rb') as f:  # Python 3: open(..., 'rb')
    dist5 = pickle.load(f)
        
with open('dists50.pkl','rb') as f:  # Python 3: open(..., 'rb')
    dist50 = pickle.load(f)

    

figure, axis = plt.subplots(dpi=600)
figure.set_size_inches(3.5, 2)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"

fsize = 10  # font size
fsize_small = 8

dist05 = [a/0.5 for a in dist05]
dist5 = [a/0.5 for a in dist5]
dist50 = [a/0.5 for a in dist50]

c3= 'red'
c1= 'gold'
c2 = 'blue'

plt.hist(dist05,bins=100,density=True,histtype='stepfilled',linewidth=0.5, alpha = 0.2, color=c1)
plt.hist(dist5,bins=100,density=True,histtype='stepfilled',linewidth=0.5, alpha = 0.2, color=c2)
plt.hist(dist50,bins=100,density=True,histtype='stepfilled',linewidth=0.5, alpha = 0.2, color=c3)

plt.hist(dist05,bins=100,density=True,histtype='step',linewidth=0.5, alpha = 0.95,  color=c1,label='$K=0.5$')
plt.hist(dist5,bins=100,density=True,histtype='step',linewidth=0.5, alpha = 0.95, color=c2,label='$K=5$')
plt.hist(dist50,bins=100,density=True,histtype='step',linewidth=0.5, alpha = 0.95, color=c3,label='$K=50$')

axis.set_xlabel('Normalized distance $l$', fontsize=fsize)
axis.set_ylabel('Density [%]', fontsize=fsize)
axis.tick_params(axis='both', which='major', labelsize=fsize_small)
axis.spines['top'].set_visible(False)
axis.spines['right'].set_visible(False)
axis.set_xlim(0,2)
#axis.set_ylim(-0.03,0.35)

print(np.var(dist05))
print(np.var(dist5))
print(np.var(dist50))

plt.tight_layout()

lgnd = plt.legend(markerscale=1,loc='upper right')

plt.savefig("histogram.pdf", format="pdf", bbox_inches="tight")