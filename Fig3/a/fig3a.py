import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
import pandas as pd 
from scipy.stats import gamma
from scipy.optimize import fsolve

################################
################################ path
################################

os.chdir(os.path.realpath(''))


################################
################################ Import Data 
################################

################################ Detection probabilities
data = pd.read_csv("../kucirka_fig2a.csv",dtype={"proba_positive": np.float32},skiprows=0) 
p_detect = data.values[:,-1]     
p_detect = np.concatenate((p_detect,np.zeros(3)))

time = np.arange(0,23+1,1)

plt.plot(time+1,p_detect,'o-',color='black',markersize=12)

plt.ylim((-0.1,1))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()






