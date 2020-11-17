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

my_data = genfromtxt('size_time.txt', delimiter=',')
#my_data = genfromtxt('size_time_poisson.txt', delimiter=',')
for i in range(len(my_data)):
    size.append(float(my_data[i][1]))
    time.append(float(my_data[i][0]))


time = np.array(time)
size = np.array(size)

################################
################################ theory
################################
########## parameters
p_hosp = 0.029

########## Hospitalization time distribution - Gamma
shape_hosp = 31.01
scale_hosp = 0.463

########## secondary infections
R0 = 2.9

########## auxiliary
dt = 0.01
tfin = 100.
t = np.arange(0,tfin+dt,dt)

################################ prediction
###### extinction probability - Poisson
def ext(x):
    return(x-np.exp(R0*(x-1)))

p_ext = fsolve(ext,0.5)
p_surv = 1 - p_ext[0]

####### epidemic size prediction
alpha = 0.210157    # Mathematica solution!
beta = 4.67873      # Mathematica solution!

######## hospitalization times
hosp_time_dist = np.zeros(len(t))

#################### computation
for k in range(200):            ### sum over geometric distribution
    
    # determine deterministic hitting time of k infected individuals
    tdet = np.log(k*alpha*beta*p_surv)/alpha
    
    if (tdet <= 0): tdet = 0
    
    # add hospitalisation time distribution to tdet
    index = int(tdet/dt)
    hosp_time_dist[index:] += geom.pmf(k,p_hosp) * gamma.pdf(t[0:(len(t)-index)],shape_hosp,scale = scale_hosp)


####### epidemic size distribution
dsize = 50
I = np.arange(dsize,15001,dsize)
dens_I = np.zeros(len(I))

ind_old = 0
for i in range(len(I)):
    tdet_index = int(np.round(np.log(I[i]*alpha*beta*p_surv)/alpha/dt))+1
        
    dens_I[i] = np.sum(hosp_time_dist[ind_old:tdet_index]*dt)/dsize
    
    ind_old = tdet_index
    
    
################################
################################ Plot
################################

plt.hist(size, bins=np.arange(0,np.max(size),50),density=True,alpha=0.5,color='C0')
plt.plot(I-dsize/2,dens_I,linewidth=5,color='C0')
plt.axvline(np.mean(size),color='C0',linewidth=3,linestyle='dotted')
plt.axvline(np.sum((I-25)*50*dens_I),color='C0',linewidth=3,linestyle='dashed')
plt.xlim((0,2000))
plt.ylim((0,0.0025))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

    
################################
################################ Size statistics
################################
print(np.mean(size))
print(np.sqrt(np.var(size)))


