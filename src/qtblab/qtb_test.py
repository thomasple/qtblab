import numpy as np
import matplotlib.pyplot as plt
from .qtb import QTB
from .atomic_units import au
from .moldyn import langevin_dynamics_1D
from .potentials import F_Mo

dt=0.2/au.fs
Tseg=1./au.ps
omegacut=9000./au.cm1
kT=300./au.kelvin
gamma0=10./au.THz
mass=au.Mprot
qtb=QTB(dt,Tseg,omegacut,kT,gamma0,mass)
print("ALL FREQUENCIES IN cm^-1, TIMES IN fs")
print(f"dt={dt*au.fs}, Tseg={Tseg*au.fs}")
print(f"domega={qtb.dom*au.cm1},first frequency {qtb.omega[1]*au.cm1}")
omegamax=np.pi/dt
n=len(qtb.omega)
print(f"omegamax={omegamax*au.cm1}, omega[-1]={qtb.omega[-1]*au.cm1}")

qtbhalf=QTB(dt,Tseg,omegacut,kT,gamma0,mass,hbar=0.5)
lgv=QTB(dt,Tseg,omegacut,kT,gamma0,mass,hbar=0)

plt.plot(qtb.omega*au.cm1,qtb.ff_kernel(),label="qtb")
plt.plot(qtb.omega*au.cm1,qtbhalf.ff_kernel(),label="qtb hbar=0.5")
plt.plot(qtb.omega*au.cm1,lgv.ff_kernel(),label="qtb hbar=0")
plt.axvline(qtb.omegacut*au.cm1,color="black")
plt.axhline(2.*mass*gamma0*kT/dt,color="red",label="classical")
plt.xlim(-2.*omegacut*au.cm1,2*omegacut*au.cm1)
plt.legend()
plt.show()

qtb.compute_spectra()
plt.plot(qtb.omega[:qtb.nom]*au.cm1,qtb.Cff)
plt.plot(qtb.omega[:qtb.nom]*au.cm1,qtb.ff_kernel()[:qtb.nom])
plt.show()

Nsteps=qtb.nseg+1
Q0=0.
P0=np.random.normal(0.,np.sqrt(mass*kT))
th=lgv
Q,P=langevin_dynamics_1D(Q0,P0,Nsteps,F_Mo,th)
print(np.sum(mass*th.Cvv)*lgv.dom*au.kelvin)
plt.plot(th.omega[:qtb.nom]*au.cm1,mass*th.Cvv,label="mCvv")
plt.plot(th.omega[:qtb.nom]*au.cm1,th.Cvf/th.gammar,label="Cvf/gamma")
plt.legend()
plt.show()
