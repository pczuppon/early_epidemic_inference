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

time = []
size = []

my_data = genfromtxt('size_time.txt', delimiter=',')
#my_data = genfromtxt('size_time_poisson.txt', delimiter=',')
for i in range(len(my_data)):
    size.append(float(my_data[i][1])-float(my_data[i][2]))
    time.append(float(my_data[i][0]))

time = np.array(time)
size = np.array(size)

np.mean(size)
np.percentile(size,5)
np.percentile(size,95)

