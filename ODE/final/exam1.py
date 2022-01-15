from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

'''
let's set
DNAdamage is y[0]
p53 is y[1]
Mdm2_cyt is y[2]
Mdm2_nuc is y[3]
'''

# all the parameters
k_prime_d53 = 0.27
k_prime2_d53 = 8.25
theta = 0.83
k_prime_d2 = 0.05
k_prime_s53 = 0.6
k_prime2_s53 = 2.56
J_s53 = 0.45
k_prime_s2 = 0.15
k_prime2_s2 = 4.23
J_s2 = 0.92
k_i = 0.41
k_o = 0.05
k_prime2_d2 = 0.79
k_repair = 0.08
J1 = 0.1
J2 = 0.1

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
    y[0] DNAdamage
    y[1] p53
    y[2] Mdm2_cyt
    y[3] Mdm2_nuc
    '''
    y0d = -k_repair*heavisible(y[0])
    k_d2 = k_prime_d2*(1+y[0])
    k_d53 = k_prime_d53+k_prime2_d53*Goldbeter_koshland(y[3],theta,J1/y[1],J2/y[1])
    y1d = k_prime_s53+k_prime2_s53*(y[2]**4/(J_s53**4+y[2]**4))-k_d53*y[1]
    y2d = k_prime_s2+k_prime2_s2*(y[1]**4/(J_s2**4+y[1]**4))-k_i*y[2]+k_o*y[3]-k_prime2_d2*y[2]
    y3d = k_i*y[2]-k_o*y[3]-k_d2*y[3]
    return [y0d,y1d,y2d,y3d]


# to reproduce column2, small DNAdamage
sol = solve_ivp(process,t_span=[0,10],y0=[0,0.19,0.19,0.78])
tmp1 = sol.y
sol = solve_ivp(process,t_span=[0,30],y0=[2.0,0.17,0.158,0.704])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.set_xlim(0,40)
ax.set_ylim(0,2.0)
ax.set_xticks([0,10,20,30,40])
ax.set_xticklabels(['40','50','60','70','80'])
ax.plot(np.arange(result.shape[1]),result[1,:],linestyle='-',color='k')
ax.plot(np.arange(result.shape[1]),result[3,:],linestyle='--',color='k')
ax1 = ax.twinx()
ax1.set_ylim(-0.5,2.5)
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey')
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam1_small_damage.pdf',bbox_inches='tight')
plt.close()

# to reproduce larger DNAdamage
sol = solve_ivp(process,t_span=[0,10],y0=[0,0.19,0.19,0.78])
tmp1 = sol.y
sol = solve_ivp(process,t_span=[0,90],y0=[5.0,0.17,0.158,0.704])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.set_xlim(0,100)
ax.set_ylim(0,2.0)
ax.set_xticks([0,20,40,60,80,100])
ax.set_xticklabels(['40','60','80','100','120','140'])
ax.plot(np.arange(result.shape[1]),result[1,:],linestyle='-',color='k')
ax.plot(np.arange(result.shape[1]),result[3,:],linestyle='--',color='k')
ax1 = ax.twinx()
ax1.set_ylim(-0.5,6)
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey')
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam1_large_damage.pdf',bbox_inches='tight')
plt.close()