from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

'''
Let's set
DNAdamage as y[0]
p53 as y[1]
Mdm2 as y[2]
'''

# parameters
k_prime_d53 = 0.005
k_prime2_d53 = 0.1
theta = 0.5
k_prime_d2 = 0.4
k_s53 = 0.0276
k_prime_s2 = 0.01
k_prime2_s2 = 0.5
J_s2 = 0.5
k_repair=0.0021429
J1 = 0.1
J2 = 0.1

#helper function
def Goldbeter_koshland(u,v,q,r):
    a = 2*u*r
    b = (v-u+v*q+u*r)**2
    c = 4*u*r*(v-u)
    d = np.sqrt(b-c)
    e = v-u+v*q+u*r+d
    return a/e

def heavisible(x):
    if x > 0:
        return 1
    else:
        return 0


def process(t,y):
    '''
    y[0] is DNAdamage
    y[1] is p53
    y[2] is Mdm2
    '''
    y0d = -k_repair*heavisible(y[0])
    k_d2 = k_prime_d2*(1+y[0])
    k_d53 = k_prime_d53+k_prime2_d53*Goldbeter_koshland(y[2],theta,J1/y[1],J2/y[1])
    y1d = k_s53-k_d53*y[1]
    y2d = k_prime_s2+k_prime2_s2*(y[2]**4/(J_s2**4+y[2]**4))-k_d2*y[2]
    return [y0d,y1d,y2d]

sol = solve_ivp(process,t_span=[0,50],y0=[0,0.36,0.71])
tmp1 = sol.y
sol = solve_ivp(process,t_span=[0,50],y0=[0.07,0.302,1.240])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.plot(np.arange(result.shape[1]),result[1,:],linestyle='-',color='k')
ax.plot(np.arange(result.shape[1]),result[2,:],linestyle='--',color='k')
ax1 = ax.twinx()
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey') # DNAdamage
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
ax.set_xticks([])
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam3_small_damage.pdf',bbox_inches='tight')
plt.close()


# larger DNA damage
sol = solve_ivp(process,t_span=[0,50],y0=[0,0.36,0.71])
tmp1 = sol.y
sol = solve_ivp(process,t_span=[0,50],y0=[0.20,0.302,1.240])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.plot(np.arange(result.shape[1]),result[1,:],linestyle='-',color='k')
ax.plot(np.arange(result.shape[1]),result[2,:],linestyle='--',color='k')
ax1 = ax.twinx()
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey') # DNAdamage
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
ax.set_xticks([])
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam3_large_damage.pdf',bbox_inches='tight')
plt.close()
