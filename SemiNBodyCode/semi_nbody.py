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
				("mu2",c_double)]


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

	def run_megno_integration(self,m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp,tFinish,dtfactor=1./30.):
		# print m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp,tFinish
		sim = Simulation()
		psim = pointer(sim)	
		self._initialize_simulation(pointer(sim),m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,ntp,ltp,etp,pomegatp)
		shortest_period = 2 * np.pi / np.max(np.array([n1,n2,ntp]))
		dt = dtfactor * shortest_period
		# print tFinish
		return self._IntegrateSimulation(psim,tFinish,dt)
		
		
#		self._free_simulation(sim)

if __name__=="__main__":
	
	w = libwrapper()
	
	m1=1.e-5
	m2=1.e-4

	e1=0.0
	e2=0.06
	pomega1=0.
	pomega2=-np.pi / 4.
	l1=0
	l2= 0.
	

	etp=0.02
	ltp=0
	pomegatp=0.



	delta1 = -0.005;
	delta2 = 0.005;
	n1 = 3. / 2. * (1+delta1)
	n2 =  2. / 3. / (1+delta2)

	tFin = 2*np.pi*3e3

	meg = w.run_megno_integration(m1,m2,n1,l1,e1,pomega1,n2,l2,e2,pomega2,1.0,ltp,etp,pomegatp,tFin,dtfactor=1./30.)
	print meg
