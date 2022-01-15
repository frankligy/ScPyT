from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

'''
Let's set
DNAdamage as y[0]
p53_active as y[1]
p53_inactive as y[2]
Mdm2_cyt as y[3]
Mdm2_nuc as y[4]
'''

# parameters
k_prime_d53 = 0.3
k_prime2_d53 = 8
theta = 0.8
k_prime_d2 = 0.1
k_prime_s53 = 0.6
k_inactivation = 0.1
k_prime_activation = 0.2
k_prime2_activation = 5
J_activation = 0.2
k_prime_s2 = 0.2
k_prime2_s2 = 3
J_s2 = 0.7
k_i = 0.4
k_o = 0.05
k_prime2_d2 = 0.7
k_repair = 0.05
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
    y[0]: DNAdamage
    y[1]: p53_active
    y[2]: p53_inactive
    y[3]: Mdm2_cyt
    y[4]: Mdm2_nuc
    '''
    p53_total = y[1]+y[2]
    k_activation = k_prime_activation+k_prime2_activation*(y[1]**3/(J_activation**3+y[1]**3))
    y0d = -k_repair*heavisible(y[0])
    k_d2 = k_prime_d2*(1+y[0])
    k_d53 = k_prime_d53+k_prime2_d53*Goldbeter_koshland(y[4],theta,J1/p53_total,J2/p53_total)
    y1d = k_activation*y[2]-k_inactivation*y[1]-k_d53*y[1]
    y2d = k_prime_s53-k_activation*y[2]+k_inactivation*y[1]-k_d53*y[2]
    y3d = k_prime_s2+k_prime2_s2*(y[1]**3/(J_s2**3+y[1]**3))-k_i*y[3]+k_o*y[4]-k_prime_d2*y[3]
    y4d = k_i*y[3]-k_o*y[4]-k_d2*y[4]
    return [y0d,y1d,y2d,y3d,y4d]

sol = solve_ivp(process,t_span=(0,50),y0=[0,0.01,0.19,0.21,0.55])
tmp1 = sol.y
sol = solve_ivp(process,t_span=(0,50),y0=[1,0.003,0.097,0.544,1.447])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.plot(np.arange(result.shape[1]),result[2,:],linestyle='-',color='k')  # p53_inactive
ax.plot(np.arange(result.shape[1]),result[4,:],linestyle='--',color='k') # Mdm2_nuc
ax1 = ax.twinx()
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey') # DNAdamage
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53_inactive','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
ax.set_xticks([])
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam2_small_damage.pdf',bbox_inches='tight')
plt.close()


# larger damage
sol = solve_ivp(process,t_span=(0,50),y0=[0,0.01,0.19,0.21,0.55])
tmp1 = sol.y
sol = solve_ivp(process,t_span=(0,50),y0=[2,0.003,0.097,0.544,1.447])
tmp2 = sol.y

result = np.concatenate([tmp1,tmp2],axis=1)

fig,ax = plt.subplots()
ax.plot(np.arange(result.shape[1]),result[2,:],linestyle='-',color='k')  # p53_inactive
ax.plot(np.arange(result.shape[1]),result[4,:],linestyle='--',color='k') # Mdm2_nuc
ax1 = ax.twinx()
ax1.plot(np.arange(result.shape[1]),result[0,:],linestyle='-',color='grey') # DNAdamage
handles = []
handles.append(mlines.Line2D([],[],linestyle='-',color='k'))
handles.append(mlines.Line2D([],[],linestyle='--',color='k'))
handles.append(mlines.Line2D([],[],linestyle='-',color='grey'))
ax.legend(handles=handles,labels=['p53_inactive','Mdm2_nuc','DNAdamage'])
ax.set_xlabel('Time',weight='bold')
ax.set_ylabel('Protein concentration',weight='bold')
ax.set_xticks([])
plt.savefig('/Users/ligk2e/Desktop/ODE/final/exam2_large_damage.pdf',bbox_inches='tight')
plt.close()




fig,ax = plt.subplots()
ax.plot(np.arange(result.shape[1]),result[1,:],linestyle='-',color='k')
ax.plot(np.arange(result.shape[1]),result[3,:],linestyle='--',color='k')
ax1 = ax.twinx()
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

