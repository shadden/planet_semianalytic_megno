from ctypes import *
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os
who =os.popen("whoami") 
if who.readline().strip() =='samuelhadden':
	SRCDIR = "/Users/samuelhadden/02_ThreePlanet_Stability/03_Semianalytic_Code/SemiNBodyCode"
else:
	SRCDIR = "/projects/p20783/sjh890/02_Chaos_Project/SemiAnalyticCodes/planet_semianalytic_megno/SemiNBodyCode"
who.close()

def get_ctype_ptr(dtype,dim,**kwargs):
	return np.ctypeslib.ndpointer(dtype=dtype,ndim=dim,flags='CONTIGUOUS',**kwargs)	
p1d=get_ctype_ptr(np.float,1)
p1dInt = get_ctype_ptr(c_int,1)


class PhaseState(Structure):
	_fields_ = [("x",c_double),
				("y",c_double),
				("vx",c_double),
				("vy",c_double),
				("dx",c_double),
				("dy",c_double),
				("dvx",c_double),
				("dxy",c_double),
				("dax",c_double),
				("day",c_double)]
class PhaseStateSimple(Structure):
	_fields_ = [("x",c_double),
				("y",c_double),
				("vx",c_double),
				("vy",c_double)
				]
class MEGNO_Auxilary_Variables(Structure):
	_fields_ = [("W",c_double),
				("Y",c_double),
				("megno",c_double)]
class Simulation(Structure):
	_fields_ = [("test_particle",POINTER(PhaseState)),
				("inner_planet",POINTER(PhaseStateSimple)),
				("outer_planet",POINTER(PhaseStateSimple)),
				("megno_aux",POINTER(MEGNO_Auxilary_Variables)),
				("mu1",c_double),
				("mu2",c_double),
				("t",c_double)]


class libwrapper(object):

	def __init__(self):
		
		self.lib = CDLL("%s/libSemiNbody.so"%SRCDIR)
		
		self._initialize_simulation=self.lib.intialize_simulation
		self._initialize_simulation.argtypes = [POINTER(Simulation)] + [c_double for i in range(4*3+2)]
		self._initialize_simulation.restype = None

		self._free_simulation=self.lib.free_simulation
		self._free_simulation.argtypes = [POINTER(Simulation)]
		self._free_simulation.restype = None

		self._free_simulation=self.lib.free_simulation
		self._free_simulation.argtypes = [POINTER(Simulation)]
		self._free_simulation.restype = None
		
		self._IntegrateSimulation = self.lib.IntegrateSimulation
		self._IntegrateSimulation.argtypes = [POINTER(Simulation), c_double, c_double]
		self._IntegrateSimulation.restype = c_double

		self._IntegrateSimulationToTime = self.lib.IntegrateSimulationToTime
		self._IntegrateSimulationToTime.argtypes = [POINTER(Simulation), c_double,c_double]
		self._IntegrateSimulationToTime.restype = None

	def run_megno_integration(self,m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp,tFinish,dtfactor=1./30.):
		sim = Simulation()
		psim = pointer(sim)	
		self._initialize_simulation(pointer(sim),m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp)
		shortest_period = 2 * np.pi / np.max(np.array([n1,n2,ntp]))
		dt = dtfactor * shortest_period
		# print tFinish
		return self._IntegrateSimulation(psim,tFinish,dt)
	def setup_integration(self,m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp,tFinish,dtfactor=1./30.):
		sim = Simulation()
		psim = pointer(sim)	
		self._initialize_simulation(psim,m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp)
		shortest_period = 2 * np.pi / np.max(np.array([n1,n2,ntp]))
		dt = dtfactor * shortest_period
		
		return psim,dt
	def Mengo_And_tLy_Integrate(self,sim,dt,tStop,Nout,MAX_MEGNO=-1):
		megnos=np.zeros(Nout)
		times = np.linspace(0,tStop,Nout)
		for i,t in enumerate(times):
			self._IntegrateSimulationToTime(pointer(sim),t,dt)
			megnos[i] =  sim.megno_aux.contents.megno
			if MAX_MEGNO > 0 and megnos[i] > MAX_MEGNO:
				break
		if i+1==Nout:
			Nfit=int(np.floor(Nout/2))
			tfit=times[Nfit:]
			megnofit=megnos[Nfit:]
			tLy= 1. / np.linalg.lstsq(np.vstack((  tfit,   np.ones(len(tfit))   )).T,megnofit)[0][0]
		else:
			tLy= 1. / np.linalg.lstsq(np.vstack((times[:i+1],np.ones(i+1))).T,megnos[:i+1])[0][0]
		return megnos[i],tLy

def get_orbital_elements(psim):
	t=psim.contents.t
	x=psim.contents.test_particle.contents.x
	x=psim.contents.test_particle.contents.x
	y=psim.contents.test_particle.contents.y
	vx=psim.contents.test_particle.contents.vx
	vy=psim.contents.test_particle.contents.vy
	
	rsq = x*x + y*y
	rinv = 1./np.sqrt(rsq)
	vsq = vx*vx + vy*vy
	angmo = x*vy - y*vx
	sma = 1.0  / (2.0 * rinv - vsq)
	ecc = np.sqrt( 1.0 - angmo * angmo / sma)
	ex = angmo * vy - x*rinv
	ey = -angmo * vx - y*rinv
	return t,sma,ecc,ex,ey
			
		
#		self._free_simulation(sim)

if __name__=="__main__":
	
	w = libwrapper()
	
	m1=1.e-5
	m2=3.e-6

	e1=0.06
	e2=0.06
	pomega1=0.
	pomega2=np.pi 
	l1=0
	l2= 0.
	

	etp=0.02
	ltp=0
	pomegatp=0.



	delta1 = -0.005;
	delta2 = 0.005;
	n1 = 3. / 2. * (1+delta1)
	n2 =  2. / 3. / (1+delta2)

	tFin = 2*np.pi*5e4
	
	psim,dt=w.setup_integration(m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,1.0,ltp,etp,pomegatp,tFin,dtfactor=1./30.)
	op=w.Mengo_And_tLy_Integrate(psim.contents,dt,tFin,10)
	print op
	meg=w.run_megno_integration(m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,1.0,ltp,etp,pomegatp,tFin,dtfactor=1./30.)
	print meg
	psim,dt=w.setup_integration(m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,1.0,ltp,etp,pomegatp,tFin,dtfactor=1./30.)
	Npts=200
	pts = np.zeros((5,Npts))
	for i,t in enumerate(np.linspace(0,tFin,Npts)):

		w._IntegrateSimulationToTime(psim,t,dt)
		pts[:,i]=get_orbital_elements(psim)
