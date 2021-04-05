import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt

'''
Figure1A linear signal response
dR/dt = k0 + k1*S - k2*R
'''

def model(R,t,S):
    k0 = 0.01
    k1 = 1
    k2 = 5
    dRdt = k0 + k1*S - k2*R
    return dRdt

R0 = [0,0.3,0.5]
t = np.linspace(0,1,10)
S = 1
result = odeint(model,R0,t,args=(S,))

# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.3')
ax.plot(t,result[:,2],label='R0=0.5')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('R')
ax.axhline(y=0.202,xmin=0,xmax=1,linestyle='--',c='k')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_solution_curve.pdf',bbox_inches='tight')
plt.close()

# rate curve
k0 = 0.01
k1 = 1
k2 = 5
fig,ax = plt.subplots()
S_options = [1,2,3]
for S in S_options:
    R = np.linspace(0,1,10)
    removal_rate = k2 * R
    production_rate = [k1 * S] * len(R)
    ax.plot(R,removal_rate,linestyle='-',c='k')
    ax.plot(R,production_rate,linestyle='--',c='k')
ax.set_xlim(0,1)
ax.set_ylim(0,6)
ax.set_xlabel('R')
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
k0 = 0.01
k1 = 1
k2 = 5
S = np.linspace(0,3,7)
R_ss = (k0 + k1*S) / k2
fig,ax = plt.subplots()
ax.plot(S,R_ss,linestyle='-',c='k')
ax.set_xlim(0,3)
ax.set_ylim(0,0.7)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1A_SR_curve.pdf',bbox_inches='tight')
plt.close()

'''
Figure1B, hyperbolic response
dRp/dt = k1*S(Rt-Rp) - k2*Rp
'''

def model(y,t,S):
    k1 = 1
    k2 = 1
    Rt = 1
    dydt = k1*S*(Rt-y) - k2*y
    return dydt

S = 1
Rp0 = [0,0.5,1]
t = np.linspace(0,1,10)
result = odeint(model,Rp0,t,args=(S,))


# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.5')
ax.plot(t,result[:,2],label='R0=1')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('Rp')
ax.axhline(y=0.5,xmin=0,xmax=1,linestyle='--',c='k')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1B_solution_curve.pdf',bbox_inches='tight')
plt.close()

# rate curve
k1 = 1
k2 = 1
Rt = 1
fig,ax = plt.subplots()
S_options = [2,4,8]
for S in S_options:
    Rp = np.linspace(0,1,10)
    removal_rate = k2 * Rp
    production_rate = k1 * S * (Rt - Rp)
    ax.plot(R,removal_rate,linestyle='-',c='k')
    ax.plot(R,production_rate,linestyle='--',c='k')
ax.set_xlim(0,1)
ax.set_ylim(0,6)
ax.set_xlabel('R')
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1B_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
k1 = 1
k2 = 1
Rt = 1
S = np.linspace(0,10,100)
Rp_ss = (S * Rt) / (k2/k1 + S)
fig,ax = plt.subplots()
ax.plot(S,Rp_ss,linestyle='-',c='k')
ax.set_xlim(0,10)
ax.set_ylim(0,1.1)
ax.set_xlabel('Signal(S)')
ax.set_ylabel('Response(R)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1B_SR_curve.pdf',bbox_inches='tight')
plt.close()

'''
sigmoidal curve
dRp/dt = (k1*S*(Rt-Rp)/(km1+Rt-Rp)) - k2*Rp/(km2+Rp)
'''

def model(Rp,t,S):
    k1 = 1
    k2 = 1
    Rt = 1
    km1 = 0.05
    km2 = 0.05
    dRpdt = (k1*S*(Rt-Rp)/(km1+Rt-Rp)) - k2*Rp/(km2+Rp)
    return dRpdt

S = 1
Rp0 = [0,0.3,1]
t = np.linspace(0,20,200)
result = odeint(model,Rp0,t,args=(S,))

# solution curve
fig,ax = plt.subplots()
ax.plot(t,result[:,0],label='R0=0')
ax.plot(t,result[:,1],label='R0=0.3')
ax.plot(t,result[:,2],label='R0=1')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('Rp')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_solution_curve.pdf',bbox_inches='tight')
plt.close()


# rate curve
k1 = 1
k2 = 1
Rt = 1
km1 = 0.05
km2 = 0.05
Rp = np.linspace(0,1,100)
fig,ax = plt.subplots()
for S in [0.25,0.5,1,1.5,2]:
    removal_rate = k2*Rp/(km2+Rp)
    production_rate = k1*S*(Rt-Rp)/(km1+Rt-Rp)
    ax.plot(Rp,removal_rate,linestyle='-',c='k')
    ax.plot(Rp,production_rate,linestyle='--',c='k')
ax.set_xlim(0,1)
ax.set_xlabel('Rp')
ax.set_ylim(0,2)
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
S_all = np.linspace(0,3,100)
def equation(Rp,S):
    k1 = 1
    k2 = 1
    Rt = 1
    km1 = 0.05
    km2 = 0.05
    return k1*S*(Rt-Rp)/(km1+Rt-Rp) - k2*Rp/(km2+Rp)

from scipy.optimize import fsolve
store = []
for S in S_all:
    Rp_ss = fsolve(equation,[1],args=(S,))[0]
    store.append(Rp_ss)

fig,ax = plt.subplots()
ax.plot(S_all,store,c='k')
ax.set_xlim(0,3)
ax.set_xlabel('Signal(S)')
ax.set_ylim(0,1.1)
ax.set_ylabel('Response(R_ss)')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1C_SR_curve.pdf',bbox_inches='tight')


'''
Figure1d perfectly adaptation to signal
dR/dt = k1*S - k2*X*R
dX/dt = k3*S - k4*X
'''

# solution curve
from gekko import GEKKO
k1,k2,k3,k4 = 2,2,1,1
m = GEKKO()
m.time = np.linspace(0,20,201)
S = np.array([0] * 40 + [1] * 40 + [2] * 40 + [3] * 40 + [4] * 41)
S_param = m.Param(value=S)
X_var = m.Var(value=0)
R_var = m.Var(value=0)
m.Equation(X_var.dt()==k3*S_param-k4*X_var)
m.Equation(R_var.dt()==k1*S_param-k2*X_var*R_var)
m.options.IMODE = 4
m.solve()

fig,ax1 = plt.subplots()
ax1.plot(m.time,R_var,c='k')
ax1.plot(m.time,X_var,c='g')
ax1.set_xlim(0,20)
ax1.set_xlabel('Time')
ax2 = ax1.twinx()
ax2.plot(m.time,S_param,c='r')
import matplotlib.lines as mlines
ax1.legend(handles=[mlines.Line2D([],[],c=i,linestyle='-') for i in 'kgr'],labels=['R','X','S'])
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1d_solution_curve.pdf',bbox_inches='tight')
plt.close()

# rate curve
S_options = [1,2,3]
fig,ax = plt.subplots()
for S in S_options:
    R = np.linspace(0,2,100)
    X_ss = k3*S/k4
    removal_rate = k2 * X_ss * R
    production_rate = [k1*S] * len(R)
    tmp_dict = {1:'r',2:'purple',3:'g'}
    ax.plot(R,removal_rate,c=tmp_dict[S],linestyle='-')
    ax.plot(R,production_rate,c=tmp_dict[S],linestyle='--')
ax.set_xlim(0,2)
ax.set_ylim(0,8)
ax.set_xlabel('R')
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1d_rate_curve.pdf',bbox_inches='tight')
plt.close()

'''
Figure1e mutual activation, one-way switch
'''
k0,k1,k2,k3,k4,J3,J4 = 0.4,0.01,1,1,0.2,0.05,0.05
def Goldbeter_Koshland(v1,v2,J1,J2):
    '''
    will return the equalibrium concertration of either phosphorated or unphosphorated form of E,
    v1,J1: incoming
    v2,J2: outcoming
    '''
    B = v2-v1+J1*v2+J2*v1
    equilibrium = 2*v1*J2/(B+np.sqrt(B**2-4*(v2-v1)*v1*J2))
    return equilibrium

def model(R,t,S):
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    dRdt = k0*EP + k1*S - k2*R
    return dRdt

t = np.linspace(0,10,100)
result = odeint(model,y0=1,t=t,args=(0,))

# solution curve, they don't ask me to do that, skip for now

# rate curve
R = np.linspace(0,0.7,100)
S_options = [0,8,16]
fig,ax = plt.subplots()
for S in S_options:
    removal_rate = k2*R
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    production_rate = k0*EP+k1*S
    ax.plot(R,removal_rate,c='k',linestyle='-')
    ax.plot(R,production_rate,c='k',linestyle='--')
ax.set_xlim(0,0.7)
ax.set_xlabel('R')
ax.set_ylim(0,0.6)
ax.set_ylabel('Rate')
plt.savefig('/Users/ligk2e/Desktop/ODE/figure1e_rate_curve.pdf',bbox_inches='tight')
plt.close()

# signal-response curve
S_options = np.linspace(0,15,100)
def equation(R,S):
    EP = Goldbeter_Koshland(k3*R,k4,J3,J4)
    return k0*EP + k1*S - k2*R
fig,ax = plt.subplots()
store_uplimb,store_downlimb = [],[]
for S in S_options:
    R_ss_uplimb = fsolve(func=equation,x0=[1],args=(S,))[0]
    R_ss_downlimb = fsolve(func=equation,x0=[0],args=(S,))[0]
    store_uplimb.append(R_ss_uplimb)
    store_downlimb.append(R_ss_downlimb)
ax.plot(S_options,store_uplimb)
ax.plot(S_options,store_downlimb)



R = np.linspace(-1,2,10000)
result = equation(R,S=5)
fig,ax = plt.subplots()
ax.plot(R,result)
ax.axhline(0,linestyle='--',c='k')

fsolve(func=equation,x0=[0.0],args=(5,))









