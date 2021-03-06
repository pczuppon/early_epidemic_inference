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

####### Infection distribution
shape_inf = 6.6
scale_inf = 0.833

########## Hospitalization time distribution - Gamma
shape_hosp = 31.01
scale_hosp = 0.463

########## secondary infections
R0 = 1.3
kappa = 0.57
p_inf = kappa/(kappa+R0)

########## auxiliary
dt = 0.01
da = dt
tfin = 100.
t = np.arange(0,tfin+dt,dt)
tfin_hosp = 80

################################ prediction
####### extinction probability - negbin
def ext(x):
    return(x-(p_inf/(1-(1-p_inf)*x))**kappa)

p_ext = fsolve(ext,0.5)[0]
p_surv = 1 - p_ext


####### epidemic size prediction
alpha = 0.0486829    # Mathematica solution!
beta = 5.28354      # Mathematica solution!


############### McKendrick-von Foerster approach
pext = [p_ext]
x = dt
for i in range(int(tfin/dt)):
    prob = kappa/(kappa + R0*(1-gamma.cdf(x,shape_inf,scale=scale_inf)))
    pext.append((prob/(1-(1-prob)*p_ext))**kappa)
    x += dt

pext = np.asarray(pext)

##### Solve pde

fx = np.zeros((int(tfin/dt),int(tfin/dt)))

tot = [1]

# initialize with 1 individual at time 0 with age 0
fx[0,0] = 1
k = 2
mu_adj = (1+pext[0])
inf_times = [0]

# infectiousness vector
mu = []
x = 0
for i in range(int(tfin/dt)):
    mu.append(R0*(gamma.pdf(x+dt/2, shape_inf, scale=scale_inf)))
    x = x + dt

mu = np.asarray(mu)

for i in range(int(tfin/dt)-1):
    # move up age
    fx[i+1,1:] = fx[i,0:-1]
    
    # add infections
    # adjusted mu
    if (tot[-1] >= k):
        prod = 1
        for j in range(len(inf_times)):
            prod = pext[i-inf_times[j]]*prod
        
        inf_times.append(i)
        mu_adj = (1- pext[0] * prod)/(1-prod)
        k += 1    
        
    fx[i+1,0] = np.sum(da*mu*mu_adj*fx[i,:])
    tot.append(np.sum(fx[i+1,:]))

tot = np.asarray(tot)

######## hospitalization times
hosp_time_dist = np.zeros(len(t))
hosp_time_dist2 = np.zeros(len(t))

for k in range(200):            ### sum over geometric distribution
    
    # determine deterministic hitting time of k infected individuals
    tdet = np.log(k*alpha*beta*p_surv)/alpha
    
    if (tdet <= 0): tdet = 0
    
    # add hospitalisation time distribution to tdet
    index = int(tdet/dt)
    index2 = np.min(np.where(tot >=k+1))
    hosp_time_dist[index:] += geom.pmf(k,p_hosp) * gamma.pdf(t[0:(len(t)-index)],shape_hosp,scale = scale_hosp)
    hosp_time_dist2[index2:] += geom.pmf(k,p_hosp) * gamma.pdf(t[0:(len(t)-index2)],shape_hosp,scale = scale_hosp)


####### epidemic size distribution
dsize = 10
I = np.arange(dsize,15001,dsize)
dens_I = np.zeros(len(I))
dens_I2 = np.zeros(len(I))

ind_old = 0
ind_old2 = 0
for i in range(len(I)):
    tdet_index = max(0,int(np.round(np.log(I[i]*alpha*beta*p_surv)/alpha/dt)))+1
    tdet_index2 = max(0,np.min(np.where(tot >=I[i])))
        
    dens_I[i] = np.sum(hosp_time_dist[ind_old:tdet_index]*dt)/dsize
    dens_I2[i] = np.sum(hosp_time_dist2[ind_old2:tdet_index2]*dt)/dsize
    
    ind_old = tdet_index
    ind_old2 = tdet_index2
    
    
################################
################################ Plot
################################

plt.hist(size, bins=np.arange(0,np.max(size),dsize),density=True,alpha=0.5,color='C1')
plt.plot(I-dsize/2,dens_I,linewidth=5,color='C1')
plt.plot(I-dsize/2,dens_I2,linewidth=5,color='black')
plt.axvline(np.mean(size),color='C1',linewidth=3,linestyle='dotted')
plt.axvline(np.sum((I-dsize/2)*dsize*dens_I),color='C1',linewidth=3,linestyle='dashed')
plt.axvline(np.sum((I-dsize/2)*dsize*dens_I2),color='black',linewidth=3,linestyle='dashed')
plt.xlim((0,400))
plt.ylim((0,0.02))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()
