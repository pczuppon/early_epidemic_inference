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
dt = 0.01
tfin = 40
t = np.arange(0,tfin+dt,dt)
inf_delay = np.zeros(len(t))

####### Weibull dist
shape_inf = 5.76
scale_inf = 1.04167

avg = int(np.round(gamma.mean(shape_inf,scale =scale_inf)))

########## secondary infections
R0 = 3
kappa = 0.58
p_inf = kappa/(kappa+R0)

####### extinction probability
def ext(x):
    return(x-(p_inf/(1-(1-p_inf)*x))**kappa)

p_ext = fsolve(ext,0.5)
p_surv = 1- p_ext[0]

t_det = np.arange(0,60,avg)
inf_det = np.zeros(len(t_det))
####### prediction
for i in range(len(t_det)):
        
    if (i > 1):
        inf_det[i] = R0 * (inf_det[i-1]-inf_det[i-2]) + inf_det[i-1]
    
    elif (i==0):
        inf_det[0] = 1
    
    elif (i==1):
        inf_det[1] = R0+1

############ Import from prediction from aux.cpp (runtime issue!)
inf_delay[0] = 1
inf_delay[1:] = genfromtxt('infectivity.txt', delimiter=',')


################################
################################ Import Data 
################################
repeats = 10000

################################ 
count = np.arange(0,repeats,1)
time = np.arange(0,tfin+1,1)
size = [ [0]*len(time)] * len(count)

for i in range(len(count)):
    my_data = genfromtxt('trajectory_tfin_%s_id_%s.txt' % (tfin,count[i]), delimiter=',')
    size[i] = my_data[0:tfin+1,1]

means = []
st_dev = []
low = []
up = []
for i in range(len(time)):
    s = []
    for j in range(len(count)):
        s.append(size[j][i])
    
    s_work = np.array(s)
    means.append(np.mean(s_work))
    st_dev.append(np.std(s_work))
    low.append(np.percentile(s_work,5))
    up.append(np.percentile(s_work,95))

plt.plot(time,means,'o')
plt.plot(t,np.cumsum(inf_delay)/p_surv,color='C0',linewidth=3)
plt.plot(t_det,inf_det/p_surv,linewidth=3)
#plt.plot(t,inf,linewidth = 3, color ='C0')
#plt.plot(time,low,color='C0',linewidth=3,alpha=0.5)/p_surv
#plt.plot(time,up,color='C0',linewidth=3,alpha=0.5)
plt.ylim((0,1000))
plt.fill_between(time,low,up,alpha=0.25)
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()






