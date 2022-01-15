
from scipy.optimize import fsolve

def steady_state(x,Stress):
    BAXm_F = BAX_T-x[0]-x[1]
    BCL2_F = BCL2_T-x[3]-x[1]
    return [-(k_f1 + k_f2 * x[2]) * x[0] + k_b * BAXm_F + k_b * x[1],
            k_asBAXmBCL2 * BAXm_F * BCL2_F - k_dsBAXmBCL2 * x[1] - k_b * x[1],
            k_prime_sBH3 + k_prime2_sBH3 * Stress - k_dBH3 * x[2] - k_asBH3BCL2 * x[2] * BCL2_F + k_dsBH3BCL2 * x[3],
            k_asBH3BCL2 * x[2] * BCL2_F - k_dsBH3BCL2 * x[3] - k_dBH3 * x[3],
            ]

root  = fsolve(func=steady_state,x0=[66.6,33.4,0,16],args=(0.5,))
BAXm_F_ss = BAX_T-root[0]-root[1]

# let's try to plot the one-paramter bifurcation (Figure4.F)
fig,ax = plt.subplots()
ax.set_xlim([0,1])
ax.set_ylim([0,100])
ax.set_xlabel('Stress',weight='bold')
ax.set_ylabel(r'$\mathbf{[BCL2]_T}$')
ax.axhline(y=80,c='grey',linestyle='--')
ax.plot([0,0.82],[46,100],c='k')
ax.plot([0,0.6],[47,100],c='k')
ax.annotate(s='SN1',xy=(0.52,80),xytext=(0.6,70),arrowprops={'arrowstyle':'->'})
ax.annotate(s='SN2',xy=(0.37,80),xytext=(0.2,90),arrowprops={'arrowstyle':'->'})
plt.savefig('/Users/ligk2e/Desktop/github/exercise_codes/ODE/final/paper2_F.pdf',bbox_inches='tight')




