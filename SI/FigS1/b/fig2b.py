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
####### Infection distribution
shape_inf = 6.6
scale_inf = 0.833

########## secondary infections
R0 = 1.3
kappa = 0.57
p_inf = kappa/(kappa+R0)

####### extinction probability (poisson)
def ext(x):
    return(x-(p_inf/(1-(1-p_inf)*x))**kappa)

p_ext = fsolve(ext,0.5)
p_surv = 1 - p_ext[0]

####### prediction
alpha = 0.0486829    # Mathematica solution!
beta = 5.28354      # Mathematica solution!

I_plot =  []
I2_plot =  []
dt = 0.01
tfin = 80
t = np.arange(0,tfin+dt,dt)

for i in range(len(t)):
    I_plot.append( np.exp(alpha*t[i]) / (alpha*beta*p_surv))
    I2_plot.append( np.exp(alpha*t[i]) / (alpha*beta))

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

time = np.arange(0,tfin+1,1)

plt.semilogy(time,means[0:len(time)],'o',color='C1')
plt.semilogy(t,I_plot,color='C1',linewidth=3)
#plt.semilogy(t,I2_plot,color='C0',linewidth=3)
#plt.plot(t_det,inf_det/p_surv,linewidth=3)
#plt.plot(t,inf,linewidth = 3, color ='C0')
#plt.plot(time,low,color='C0',linewidth=3,alpha=0.5)/p_surv
#plt.plot(time,up,color='C0',linewidth=3,alpha=0.5)
plt.ylim((0,5*1000))
plt.xlim((0,80))
plt.fill_between(time,low[0:len(time)],up[0:len(time)],alpha=0.25,color='C1')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()


