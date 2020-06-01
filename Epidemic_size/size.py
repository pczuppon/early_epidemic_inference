import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
from scipy.stats import gamma
from scipy.stats import geom
from scipy.optimize import fsolve

################################
################################ path
################################

os.chdir(os.path.realpath(''))

################################
################################ Import Data 
################################
time = []
size = []
init = []

#my_data = genfromtxt('size_time.txt', delimiter=',')
my_data = genfromtxt('size_time_poisson.txt', delimiter=',')
for i in range(len(my_data)):
    size.append(float(my_data[i][1]))
    time.append(float(my_data[i][0]))
    init.append(float(my_data[i][2]))

time = np.array(time)
size = np.array(size)
init = np.array(init)

################################
################################ theory
################################
########## parameters
p_hosp = 0.026

########## Hospitalization time distribution - Gamma
shape_hosp = 8.
scale_hosp = 1.25

########## secondary infections
R0 = 3.41
kappa = 0.91
p_inf = kappa/(kappa+R0)

########## auxiliary
dt = 0.01
tfin = 100.
t = np.arange(0,tfin+dt,dt)
tfin_hosp = 80

################################ prediction
####### extinction probability
def ext(x):
    return(x-(p_inf/(1-(1-p_inf)*x))**kappa)

p_ext = fsolve(ext,0.5)
p_surv = 1 - p_ext[0]

###### Poisson offspring number
def ext_p(x):
    return(x-np.exp(R0*(x-1)))

p_ext_p = fsolve(ext_p,0.5)
p_surv_p = 1 - p_ext_p[0]

######## simulated popsizes
inf_delay = np.zeros(10001)
inf_delay[0] = 1
inf_delay[1:] = genfromtxt('infectivity_1.txt', delimiter=',')

inf = np.cumsum(inf_delay)

######## hospitalization times
#t_det = genfromtxt('hosp_times.txt')
t_det = genfromtxt('hosp_times_poisson.txt')


############### size at first event
step_size = 50
dx_det = np.arange(0,3100,step_size)
dist_det = np.zeros(len(dx_det))

k = 1
last = 0
for i in range(len(t_det)):
    #if (inf[i]/p_surv > k*step_size):
    if (inf[i]/p_surv_p > k*step_size):
        #if (inf[i]/p_surv > dx_det[-2]):
        print(i)
        if (inf[i]/p_surv_p > dx_det[-2]):
            break
        else:
            dist_det[k-1] += np.sum(t_det[last:i])
            k += 1
            last = i     


################################
################################ Plot
################################

################### pop sizes
plt.hist(size, bins=np.arange(0,np.max(size),50),density=True,alpha = 0.5)
#plt.hist(size2, bins=np.arange(0,np.max(size2),50),density=True,alpha = 0.5,color='C1')
plt.plot(dx_det+step_size/2,dist_det/step_size,linewidth=5,color='C0')
#plt.plot(dx_det+step_size/2,dist_det2/step_size,linewidth=5,color='C1')
#plt.axvline(np.percentile(size,2.5),color = 'C0',linewidth=3)
#plt.axvline(np.median(size),color = 'C0',linewidth=2)
#plt.axvline(np.percentile(size,97.5),color = 'C0',linewidth=3)
#plt.axvline(np.median(size2),color='C1',linewidth=3)
#plt.axvline(276,color='C1',linewidth=3,linestyle='dashed')
plt.xlim((0,1500))
plt.ylim((0,0.004))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()


