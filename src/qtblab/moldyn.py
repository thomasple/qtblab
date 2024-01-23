import numpy as np
from numpy import pi
from numba import njit
from numpy.random import normal
import os
import sys
from copy import copy

from .atomic_units import au

def langevin_dynamics_1D(Q0,P0,Nsteps,Force,thermostat):
  Q=Q0
  P=P0
  dt=thermostat.dt
  mass=thermostat.mass
  F=Force(Q)
  temp=0.
  for istep in range(Nsteps-1):
    P+=0.5*dt*F
    Q+=0.5*dt*P/mass
    thermostat.apply_O(P)
    temp+=P**2/mass
    print(temp/(istep+1)*au.kelvin)
    Q+=0.5*dt*P/mass
    F=Force(Q)
    P+=0.5*dt*F
  return Q,P 
