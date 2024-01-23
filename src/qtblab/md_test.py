import numpy as np
import matplotlib.pyplot as plt
from numba import njit

from .qtb import QTB,adQTB
from .atomic_units import au
from .moldyn import langevin_dynamics_1D

ndof=2
mass=np.array([au.Mprot for i in range(ndof)])
w=np.array([1600.,3500.])/au.cm1
kHO=mass*w**2
c3=1.e-2
@njit
def F_CO(Q):
  return -kHO*Q[:] - 3*c3*(Q[0]-Q[1])**2 

dt=0.2/au.fs
Tseg=0.5/au.ps
omegacut=9000./au.cm1
kT=300./au.kelvin
gamma0=5./au.THz
qtbs=[QTB(dt,Tseg,omegacut,kT,gamma0,mass[i],hbar=0.) for i in range(ndof)]
adapter=adQTB(qtbs,0.0/au.cm1)
nseg=qtbs[0].nseg
Nblocks=2000
Ntherm=10
plot_stride=10
temp=0.
Q=np.zeros(ndof)
P=np.random.normal(0.,np.sqrt(mass*kT))
Force=F_CO
F=Force(Q)
temp_avg=0
for iblock in range(Nblocks):
  temp=0.
  for istep in range(nseg):
    P+=0.5*dt*F
    Q+=0.5*dt*P/mass
    P=np.array([qtbs[i].apply_O(P[i]) for i in range(ndof)])
    temp+=np.mean(P**2/mass)
    Q+=0.5*dt*P/mass
    F=Force(Q)
    P+=0.5*dt*F
  temp_avg+=temp/nseg
  print(temp/nseg*au.kelvin,temp_avg/(iblock+1)*au.kelvin)
  if(iblock>=Ntherm):
    adapter.adapt_gammar()
    if (iblock+1-Ntherm)%plot_stride==0:
      plt.clf()
      #plt.plot(adapter.omega*au.cm1,adapter.gammar*au.THz,label='gammar (adapter)')
      plt.plot(adapter.omega*au.cm1,adapter.Cff_avg/qtbs[0].ff_kernel()[:adapter.nom],label='Cff')
      #plt.plot(adapter.omega*au.cm1,adapter.mgCvv_avg,label='mgCvv')
      #plt.plot(adapter.omega*au.cm1,adapter.Cvf_avg,label='Cvf')
      plt.legend()
      plt.pause(0.01)

plt.show()
