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

################################
################################ theory
################################
dt = 0.1
tfin = 40
t = np.arange(0,tfin+dt,dt)

################################
################################ Import Data 
################################
repeats = 10000

################################ 
count = np.arange(0,repeats,1)
time = np.arange(0,tfin+1,dt)
size = [ [0]*len(time)] * len(count)

for i in range(len(count)):
    my_data = genfromtxt('trajectory_tfin_%s_id_%s.txt' % (tfin,count[i]), delimiter=',')
    size[i] = my_data[0:len(time)+1,1]

means = []
for i in range(len(time)):
    s = []
    for j in range(len(count)):
        s.append(size[j][i])
    
    s_work = np.array(s)
    means.append(np.mean(s_work))

plt.plot(time,means,'o')
#plt.plot(time,low,color='C0',linewidth=3,alpha=0.5)/p_surv
#plt.plot(time,up,color='C0',linewidth=3,alpha=0.5)
plt.ylim((0,1000))
#plt.fill_between(time,low,up,alpha=0.25)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()


data = open('means.txt','w')
for i in range(len(means)):
    data.write(str(means[i]))
    data.write("\n")

data.close()



