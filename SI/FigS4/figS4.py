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
repeats = 10000
tfin = 80

################################ 
count = np.arange(0,repeats,1)
time = np.arange(0,tfin+1,1)
size = [ [0]*len(time)] * len(count)

for i in range(len(count)):
    my_data = genfromtxt('trajectory_tfin_%s_id_%s.txt' % (tfin,count[i]), delimiter=',')
    size[i] = my_data[0:tfin+1,1]

means = []
low = []
up = []
for i in range(len(time)):
    s = []
    for j in range(len(count)):
        s.append(size[j][i])
    
    s_work = np.array(s)
    means.append(np.mean(s_work))
    low.append(np.percentile(s_work,5))
    up.append(np.percentile(s_work,95))

plt.plot(time,means[0:len(time)],'o',color='C0')
plt.plot(t,I_plot,color='C0',linewidth=3)
plt.fill_between(time,low[0:len(time)],up[0:len(time)],alpha=0.25,color='C0')
################################
################################ theory
################################
########## parameters
p_hosp = 0.029

########## Hospitalization time distribution - Gamma
shape_hosp = 31.
scale_hosp = 0.463

########## secondary infections
R0 = 1.3

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
alpha = 0.0486829   # Mathematica solution!
beta = 5.28354     # Mathematica solution!

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
I = np.arange(dsize,2001,dsize)
dens_I = np.zeros(len(I))

ind_old = 0
for i in range(len(I)):
    tdet_index = int(np.round(np.log(I[i]*alpha*beta*p_surv)/alpha/dt))+1
        
    dens_I[i] = np.sum(hosp_time_dist[ind_old:tdet_index]*dt)/dsize
    
    ind_old = tdet_index

################ initial cases (average)  
I0 = 1/p_hosp

################ we  
hosp_time_avg = np.sum(hosp_time_dist*t*dt)
I_we = np.sum((I-25)*50*dens_I)

################################
################################ Import size and time data 
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

hosp_time_avg = np.mean(time)
################ rule of thumb (Jombart)
ps_I0 = 1- p_ext[0]**I0
Ihosp = I0*np.exp(alpha*shape_hosp*scale_hosp)/(alpha*beta*ps_I0)

################################
################################ Plot
################################
tplot = np.arange(0,hosp_time_avg+40,dt)
avg_hosp = scale_hosp*shape_hosp

tplot_Jombart = np.arange(hosp_time_avg-avg_hosp,hosp_time_avg+40,dt)

plt.plot(tplot,np.exp(tplot*alpha)/(alpha*beta*p_surv),linewidth=3,color='C0')
plt.plot(tplot_Jombart,I0*np.exp(alpha*tplot[0:len(tplot_Jombart)])/(alpha*beta*ps_I0),color='black',linewidth=3)
#plt.axvline(hosp_time_avg,color='C0',linewidth=2,linestyle='dashed')
plt.axvline(np.mean(time),color='C0',linewidth=2,linestyle='dashed')
#plt.plot(t,[np.mean(size)]*len(t),color='C0',linewidth=2,linestyle='dotted')
plt.plot(t,[1/p_hosp]*len(t),color='black',linewidth=2,linestyle='dotted')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.ylim((0,500))
plt.xlim((0,80))
plt.show()

