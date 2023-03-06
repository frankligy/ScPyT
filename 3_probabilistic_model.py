import numpy as np
from scipy.stats import norm,uniform,invgamma
import seaborn as sns
import matplotlib.pyplot as plt

norm.pdf(2,loc=0,scale=1)
norm.cdf(0,loc=0,scale=1)
norm.ppf(0.5,loc=0,scale=1)
norm.sf(0,loc=0,scale=1)
norm.isf(0.5,loc=0,scale=1)
norm.rvs(size=100,loc=0,scale=1,random_state=1)

# accept-reject method
x = np.linspace(0,1,1000)
f_x = x**1.7*((1-x)/(1+x))**5.3
M = f_x.max()
def f(x):
    return x**1.7*((1-x)/(1+x))**5.3
def g(x):
    return uniform.pdf(x,loc=0,scale=1)

n = 2500
u = uniform.rvs(loc=0,scale=1,size=n,random_state=1)
y = uniform.rvs(loc=0,scale=1,size=n,random_state=2)
f_y = f(y)
g_y = g(y)
acc_or_rej = u <= f_y / (M * g_y)
accepted_y = y[acc_or_rej]
sns.histplot(accepted_y)

# MCMC
def f(x):
    import math
    return 2*x**2*(1-x)**8*math.cos(4*math.pi*x)**2
def q(x,y):
    return norm.pdf(x,loc=y,scale=0.1)

n = 10000
x = np.zeros(n)
x[0] = norm.rvs(loc=0,scale=0.03,size=1)[0]
for i in range(n-1):
    while True:
        x_cand = norm.rvs(loc=x[i],scale=0.03,size=1)[0]
        if x_cand >= 0 and x_cand <= 1:
            break
    if x_cand >= 0 and x_cand <= 1:
        rho = (q(x[i],x_cand)/q(x_cand,x[i]))*(f(x_cand)/f(x[i]))
        alpha = min(1,rho)
        u = uniform.rvs(loc=0,scale=1,size=1)[0]
        if u < alpha:
            x[i+1] = x_cand
        else:
            x[i+1] = x[i]
sns.histplot(x)
fig,ax = plt.subplots()
ax.plot(np.arange(10000),x)

# Gibbs sampling and Hierarchical bayesian
nruns = 10000
K = 100 # K = 100, K hospitals
n = 1000 # each hospital has 1000 patients
y = np.zeros(K,n)  # generated observations
lambda_est = np.zeros(K,nruns)
sigma_est = np.zeros(nruns)
mu_est = np.zeros(nruns)
tau_est = np.zeros(nruns)

for i in range(K):
    loc = uniform.rvs(loc=0,scale=10,size=1)[0]
    scale = uniform.rvs(loc=0,scale=0.1,size=1)[0]
    lambda_est[i,0] = norm.rvs(loc=loc,scale=scale,size=1)[0]
sigma_est[0] = uniform.rvs(loc=0,scale=0.1,size=1)[0]
mu_est[0] = norm.rvs(loc=uniform.rvs(loc=0,scale=10,size=1)[0],scale=1,size=1)[0]
tau_est[0] = uniform.rvs(loc=0,scale=0.1,size=1)[0]

for runs in range(1,n_runs-1,1):
    # estimate lambda
    for i in range(1,K-1):
        mean = i/math.sqrt(1/tau_est[runs-1]) + n/sigma_est[runs-1]
        std = (mean^2)*(mu_est[runs-1]/(tau_est[runs-1])+y[i,:].mean()*n/sigma_est[runs-1])
        lambda_est[i,runs] = norm.rvs(loc=mean,scale=std,size=1)[0]
    # estimate sigma
    sigma_sum_term = 0
    for i in range(K):
        for j in range(n):
            sigma_sum_term += (y[i,j]-lambda_est[i,runs])**2
    sigma_est[runs] = invgamma(loc=K*n/2,scale=sigma_sum_term/2)
    # estimate tau
    tau_sum_term = 0
    for i in range(K):
        tau_sum_term += (lambda_est[i,runs]-mu_est[runs-1])**2
    tau_est[runs] = invgamma(loc=K/2,scale=tau_sum_term/2)
    # estimate mu
    mu_est[runs] = norm.rvs(loc=lambda_est[:,runs-1].mean(),scale=math.sqrt(tau_est[runs]/2))


