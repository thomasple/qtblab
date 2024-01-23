import numpy as np

from .atomic_units import au
from .counter import Counter


class QTB:
  def __init__(self,dt,Tseg,omegacut,kT,gamma0,mass,hbar=1.):
    self.dt=dt
    self.kT=kT
    self.omegacut=omegacut
    self.hbar = hbar
    self.mass=mass

    self.omegasmear = omegacut*0.04 
    self.nseg=int(Tseg/dt)
    self.Tseg=self.nseg*dt
    self.dom= 2.*np.pi/(3.*self.Tseg)
    self.nom=int(omegacut/self.dom)
    self.omega = self.dom*np.arange((3*self.nseg)//2+1)

    if omegacut>=self.omega[-1] :
      raise ValueError("omegacut must be smaller than pi/dt")
    self.u = 0.5*hbar*np.abs(self.omega)/kT
    self.theta=kT*np.ones_like(self.omega)
    if hbar>0:
      self.theta[1:]*=self.u[1:]/np.tanh(self.u[1:])
    self.cutoff = 1./(1.+np.exp((np.abs(self.omega)-self.omegacut)/self.omegasmear))
    

    if gamma0>0.5*omegacut:
      raise ValueError("gamma0 must be much smaller than omegacut (at least half)")
    self.gamma0=gamma0
    self.gammar=gamma0*np.ones(self.nom)
    self.a1=np.exp(-gamma0*dt)
    self.OUcorr = ( ( 1. - 2.*self.a1*np.cos(self.omega*dt) + self.a1**2 )
                      / (dt**2*(gamma0**2+self.omega**2)) 
                  )
    self.white_noise = np.random.normal(0.,1.,(3*self.nseg))
    self.refresh_force()
    self.vel=np.zeros((self.nseg))
    self.counter = Counter(self.nseg)

  def ff_kernel(self):
    kernel = self.theta*self.cutoff*self.OUcorr
    kernel[:self.nom]*=self.gammar/self.gamma0
    return kernel*2*self.mass*self.gamma0/self.dt
  
  def refresh_force(self):
    self.white_noise[:2*self.nseg]=np.copy(self.white_noise[self.nseg:])
    self.white_noise[2*self.nseg:]=np.random.normal(0.,1.,self.nseg)
    s=np.fft.rfft(self.white_noise)*np.sqrt(self.ff_kernel())
    self.force = np.fft.irfft(s,3*self.nseg)[self.nseg:2*self.nseg]
  
  def apply_O(self,p):
    istep=self.counter.i()
    force = self.force[istep]
    self.vel[istep] = (p*np.sqrt(self.a1) + 0.5*self.dt*force)/self.mass
    self.counter.increment()
    if self.counter.is_reset_step():
      self.compute_spectra() 
      self.refresh_force()
    return p*self.a1 + self.dt*force
    
    
  def compute_spectra(self):
    sf=np.fft.rfft(self.force,3*self.nseg,norm="ortho")
    sv=np.fft.rfft(self.vel,3*self.nseg,norm="ortho")
    self.Cvv = (np.abs(sv)**2)[:self.nom]*self.dt
    self.Cvf = np.real(sv*np.conj(sf))[:self.nom]*self.dt
    self.Cff = (np.abs(sf)**2)[:self.nom]*self.dt

class adQTB:
  def __init__(self,qtbs,Agamma,gammar=None):
    self.qtbs=qtbs
    self.nom=qtbs[0].nom
    self.omega=qtbs[0].omega[:self.nom]
    self.ndof=len(qtbs)
    self.Agamma = Agamma
    self.a1 = Agamma*qtbs[0].Tseg*qtbs[0].gamma0
    self.gammar = np.copy(qtbs[0].gammar)
    if gammar is not None:
      self.gammar[:]=np.copy(gammar)
    for qtb in qtbs:
      qtb.gammar = self.gammar
    self.dFDT = np.random.normal(0.,self.Agamma,self.nom)
    self.dFDT_avg = np.zeros_like(self.dFDT)
    self.mgCvv_avg = np.zeros_like(self.dFDT)
    self.Cvf_avg = np.zeros_like(self.dFDT)
    self.Cff_avg = np.zeros_like(self.dFDT)
    self.nadapt = 0
  
  def compute_dFDT(self):
    self.mgCvv = np.zeros_like(self.qtbs[0].Cvv)
    self.Cvf = np.zeros_like(self.qtbs[0].Cvf)
    self.Cff = np.zeros_like(self.qtbs[0].Cff)
    for qtb in self.qtbs:
      self.mgCvv += qtb.mass*qtb.Cvv*self.gammar/self.ndof
      self.Cvf += qtb.Cvf/self.ndof
      self.Cff += qtb.Cff/self.ndof
    self.dFDT = self.mgCvv-self.Cvf

  def update_gammar(self):
    self.gammar-= self.a1*self.dFDT/np.linalg.norm(self.dFDT)
    for qtb in self.qtbs:
      qtb.gammar = self.gammar
  
  def adapt_gammar(self):
    self.nadapt+=1
    self.compute_dFDT()
    self.dFDT_avg += (self.dFDT-self.dFDT_avg)/self.nadapt
    self.mgCvv_avg += (self.mgCvv-self.mgCvv_avg)/self.nadapt
    self.Cvf_avg += (self.Cvf-self.Cvf_avg)/self.nadapt
    self.Cff_avg += (self.Cff-self.Cff_avg)/self.nadapt
    self.update_gammar()
