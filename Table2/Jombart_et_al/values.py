import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
from scipy.stats import gamma
from scipy.optimize import fsolve

################################
################################ path
################################

os.chdir(os.path.realpath(''))

size = []

my_data = genfromtxt('size.txt', delimiter=',')
for i in range(len(my_data)):
    size.append(float(my_data[i]))

size = np.array(size)

np.mean(size)
np.percentile(size,5)
np.percentile(size,95)

