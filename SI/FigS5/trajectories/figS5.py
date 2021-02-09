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
R0 = 1.1

####### extinction probability (poisson)
def ext(x):
    return(x-np.exp(R0*(x-1)))

p_ext = fsolve(ext,0.5)[0]
p_surv = 1 - p_ext

####### prediction
alpha = 0.01746    # Mathematica solution!
beta = 5.41904      # Mathematica solution!
I_plot =  []
dt = 0.01
da = 0.01
tfin = 100
t = np.arange(0,tfin+dt,dt)

for i in range(len(t)):
    I_plot.append( np.exp(alpha*t[i]) / (alpha*beta*p_surv))


############### McKendrick-von Foerster approach
###### age dependent survival probability
pext = [p_ext]
x = dt
for i in range(int(tfin/dt)):
    pext.append(p_ext*np.exp((1-p_ext)*R0*gamma.cdf(x,shape_inf,scale=scale_inf)))
    x += dt

pext = np.asarray(pext)


##### PDE solution
fx = np.zeros((int(tfin/dt)+1,int(tfin/dt)+1))

tot = [1]

# initialize with 1 individual at time 0 with age 0
fx[0,0] = 1
k = 2
mu_adj = (1+pext[0])
inf_times = [0]

# infectiousness vector
mu = []
x = 0
for i in range(int(tfin/dt)+1):
    mu.append(R0*(gamma.pdf(x+dt/2, shape_inf, scale=scale_inf)))
    x = x + dt

mu = np.asarray(mu)

for i in range(int(tfin/dt)):
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

plt.semilogy(time,means[0:len(time)],'o',color='C0')
plt.semilogy(t,I_plot,color='C0',linewidth=3)
plt.plot(t,tot,color='black',linewidth=3)
plt.ylim((0,1000))
plt.xlim((0,tfin))
plt.fill_between(time,low[0:len(time)],up[0:len(time)],alpha=0.25,color='C0')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()






