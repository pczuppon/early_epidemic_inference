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

means = []
data = genfromtxt('means.txt')
for i in range(len(data)):
    means.append(data[i])

######## hospitalization times
t_det = np.zeros(len(t))

s= 1
j = 0
while (s <= 250):
    print(s)
    
    ############# find time where to start with s individuals initially 
    #while (inf[j]/p_surv < s):
    while (inf[j]/p_surv_p < s):
        j += 1
    
    #while (means[j] < s):
    #    j += 1
    
    ############# update time proba 
    k = 1
    while (k <= tfin_hosp/dt):
        if (j + k < len(t)):
            t_det[j+k] += (gamma.cdf(t[k],shape_hosp,scale=scale_hosp)-gamma.cdf(t[k-1],shape_hosp,scale=scale_hosp)) * geom.pmf(s,p_hosp)
            k += 1
        else:
            k = tfin_hosp/dt + 1
    s += 1


################################
################################ Plot
################################
##### popsize at detection
plt.hist(init, bins=np.arange(0,np.max(init),2),density=True,alpha=0.5)
nplot = np.arange(1,1000,1)
plt.plot(nplot,geom.pmf(nplot,p_hosp))
plt.xlim((0,400))
#plt.ylim((0,0.001))
plt.show()

################## hospitalization times
### fitting a gamma
shapefit, locfit, scalefit = gamma.fit(time,floc=0)

plt.hist(time, bins=np.arange(0,np.max(time),2),density=True,alpha=0.5)
plt.plot(t,t_det/dt,linewidth=3,color='C0')
#plt.plot(t,gamma.pdf(t,shapefit,scale=scalefit),color='C1',linewidth=3)
#plt.plot(t,gamma.pdf(t,np.round(shapefit),scale=scalefit),color='C2',linewidth=3)
plt.xlim((0,tfin_hosp))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

###### write hosp-time distribution in file
#data = open('hosp_times_poisson.txt','w')
#for i in range(len(t_det)):
#    data.write(str(t_det[i]))
#    data.write("\n")
#
#data.close()
