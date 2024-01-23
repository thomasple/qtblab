from math import exp,sqrt
from atomic_units import au
import numpy as np
import numba

D=20./au.kcalperMol
alpha=2.5*au.bohr
xmax=2.5/au.bohr
QO=0.1
asym=0.707
d0=0.95/au.bohr
ADM=2.32e5/au.kcalperMol
BDM=3.15*au.bohr
CDM=2.31e4/(au.kcalperMol*au.bohr**6)

@numba.jit(nopython=True)
def Hessian_Mo(x):	
	if x<=xmax:
		return D*2.*alpha*alpha*(2.*exp(-2.*alpha*x)
				-exp(-alpha*x))
	else:
		return D*2.*alpha*alpha*(2*exp(-2.*alpha*x) 
				-exp(-alpha*x))+D*1.*exp(1.*(x-xmax))

@numba.jit(nopython=True)
def F_Mo(x):	
	if x<=xmax:
		return -D*(-2.*alpha*exp(-2.*alpha*x)
				+2.*alpha*exp(-alpha*x))
	else:
		return -D*(-2.*alpha*exp(-2.*alpha*x) 
				+2.*alpha*exp(-alpha*x)+1.*exp(1.*(x-xmax)))

@numba.jit(nopython=True)
def F_QO(x):	
	return -4*QO*x**3

@numba.jit(nopython=True)
def F_OHO(X):
	if X.shape != (2,):
		return None
	F=np.zeros(2)
	F[1]=-D*2*alpha*( 
				-exp(-2*alpha*(X[2]/2+X[1]-d0)) 
				+ exp(-alpha*(X[2]/2+X[1]-d0)) 
				+ asym*( 
					exp(-2*alpha*(X[2]/2-X[1]-d0)/asym) 
					- exp(-alpha*(X[2]/2-X[1]-d0)/asym) 
				) )
	F[2]=-D*alpha*( 
				-exp(-2*alpha*(X[2]/2+X[1]-d0)) 
				+exp(-alpha*(X[2]/2+X[1]-d0)) 
				+asym*( 
					-exp(-2*alpha*(X[2]/2-X[1]-d0)/asym) 
					+exp(-alpha*(X[2]/2-X[1]-d0)/asym) 
				)) - ADM*BDM*exp(-BDM*X[2])+6*CDM/(X[2]**7)

	return F
	
			