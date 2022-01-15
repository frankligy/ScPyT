from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

# parameter
BCL2_T = 80
BAX_T = 100
Stress = 0.6
k_f1 = 1
k_f2 = 3
k_b = 2
k_asBAXmBCL2 = 90
k_dsBAXmBCL2 = 0.05
k_asBH3BCL2 = 10
k_dsBH3BCL2 = 0.01
k_prime_sBH3 = 0.1
k_prime2_sBH3 = 0.6
k_dBH3 = 0.01

# process
def process(t,y):
    '''
    y[0] is BAX
    y[1] is BAXm/BCL2
    y[2] is BH3_F
    y[3] is BH3/BCL2
    '''
    BAXm_F = BAX_T-y[0]-y[1]
    BCL2_F = BCL2_T-y[3]-y[1]
    y0d = -(k_f1 + k_f2 * y[2]) * y[0] + k_b * BAXm_F + k_b * y[1]
    y1d = k_asBAXmBCL2 * BAXm_F * BCL2_F - k_dsBAXmBCL2 * y[1] - k_b * y[1]
    y2d = k_prime_sBH3+k_prime2_sBH3*Stress-k_dBH3*y[2]-k_asBH3BCL2*y[2]*BCL2_F+k_dsBH3BCL2*y[3]
    y3d = k_asBH3BCL2*y[2]*BCL2_F-k_dsBH3BCL2*y[3]-k_dBH3*y[3]
    return [y0d,y1d,y2d,y3d]

sol = solve_ivp(process,t_span=(0,400),y0=[66.6,33.4,0,16])
tmp = sol.y

BH3_T = np.sum(tmp[2:,:],axis=0)
BAXm_T = BAX_T - tmp[0,:]
BAXm_F = BAX_T - np.sum(tmp[0:2,:],axis=0)
BH3_F = tmp[2,:]

fig,ax = plt.subplots()
x = np.arange(tmp.shape[1])
ax.plot(x,BH3_T,linestyle='--',color='k',label=r'$[BH3]_{T}$')
ax.plot(x,BAXm_T,linestyle='-',color='k',label=r'$[BAXm]_{T}$')
ax.plot(x,BAXm_F,linestyle='-',color='grey',label=r'$[BAXm]_{F}$')
ax.plot(x,BH3_F,linestyle='--',color='grey',label=r'$[BH3]_{F}$')
ax.legend()
ax.set_xticks([])
ax.set_xlabel('t,min')
ax.set_ylabel('Levels')
plt.savefig('/Users/ligk2e/Desktop/ODE/final/paper2.pdf',bbox_inches='tight')
plt.close()

