from ctypes import *
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os
who =os.popen("whoami") 
if who.readline().strip() =='samuelhadden':
	SRCDIR = "/Users/samuelhadden/02_ThreePlanet_Stability/03_Semianalytic_Code"
else:
	SRCDIR = "/projects/p20783/sjh890/02_Chaos_Project/SemiAnalyticCodes/planet_semianalytic_megno"
who.close()

def get_ctype_ptr(dtype,dim,**kwargs):
	return np.ctypeslib.ndpointer(dtype=dtype,ndim=dim,flags='CONTIGUOUS',**kwargs)	
p1d=get_ctype_ptr(np.float,1)
p1dInt = get_ctype_ptr(c_int,1)

# typedef struct ResonanceData {
# 	int Nres,MaxOrder;
# 	int* ResonanceIndices;
# 	double* ResonanceCoefficients;
# } ResonanceData;

# 	double alphaIn,alphaOut;
# 	double fSecIn,gSecIn;
# 	double fSecOut,gSecOut;

def get_fcoeff_ratio(pratio):
	return np.power(pratio , -0.552172)

# Resonant terms to include
def full_resonance_list(resJ,order):
	return [[resJ,order,i] for i in range(0,order+1)]

def Farey_Sequence_N(n):
	# Farey sequence function lifted from Wikipedia
	sequence = []
	a,b,c,d=0,1,1,n
	sequence.append((a,b))
	while c <= n:
		k = int( (n+b) / d )
		a,b,c,d = c , d, (k*c - a) , (k*d-b)
		sequence.append((a,b))
	return sequence
	

class ResonanceData(Structure):
	_fields_ = [("Nres",c_int),("IncludeZeroth",c_int),
				("MaxOrder",c_int),
				("ResonanceIndices",POINTER(c_int)),
				("ResonanceCoefficients",POINTER(c_double))]
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
class ActionAnglePhaseState(Structure):
	_fields_ = [("L",c_double),
				("l",c_double),
				("Y",c_double),
				("X",c_double),
				("dL",c_double),
				("dl",c_double),
				("dY",c_double),
				("dX",c_double),
				("dLdot",c_double),
				("dldot",c_double),
				("dYdot",c_double),
				("dXdot",c_double)]
class SimulationParameters(Structure):
	_fields_ = [("mu1",c_double),
				("mu2",c_double),
				("n1",c_double),
				("n2",c_double),
				("e1",c_double),
				("e2",c_double),
				("lambda2",c_double),
				("varpi2",c_double),
				("alphaIn",c_double),
				("alphaOut",c_double),
				("fSecIn",c_double),
				("gSecIn",c_double),
				("fSecOut",c_double),
				("gSecOut",c_double)
				]
class MEGNO_Auxilary_Variables(Structure):
	_fields_ = [("W",c_double),
				("Y",c_double),
				("megno",c_double)]
class ActionAngleSimulation(Structure):
	_fields_ = [		("state", ActionAnglePhaseState),
				("parameters",SimulationParameters),
				("megno",MEGNO_Auxilary_Variables),
				("rIn",ResonanceData),
				("rOut",ResonanceData),
				("t",c_double),
				("dt",c_double)
				]


class libwrapper(object):

	def __init__(self):
		
		self.lib = CDLL("%s/libsemianalyticMEGNO.so"%SRCDIR)
		self._MEGNO_Integration = self.lib.MEGNO_Integration
		self._MEGNO_Integration.argtypes = [c_double for i in range(7)]
		self._MEGNO_Integration.restype = c_double
		
		self._CircularFirstOrderResonanceMEGNOIntegration = self.lib.CircularFirstOrderResonanceMEGNOIntegration
		self._CircularFirstOrderResonanceMEGNOIntegration.argtypes = [c_int,p1dInt,c_int,p1dInt,c_double,c_double,c_double,c_double,c_double]
		self._CircularFirstOrderResonanceMEGNOIntegration.restype = c_double
		
		self._ActionAnglePhaseStateInitialize = self.lib.ActionAnglePhaseStateInitialize
		self._ActionAnglePhaseStateInitialize.argtypes = [POINTER(ActionAnglePhaseState),c_double,c_double,c_double,c_double]
		self._ActionAnglePhaseStateInitialize.restype = None

		self._InitializeActionAngleSimulation = self.lib.InitializeActionAngleSimulation
		self._InitializeActionAngleSimulation.argtypes = [POINTER(ActionAngleSimulation),c_int,c_int,p1dInt,c_int,c_int,p1dInt] + [c_double for i in range(13)]
		self._InitializeActionAngleSimulation.restype = None

		self._SimulationStep = self.lib.SimulationStep
		self._SimulationStep.argtypes = [POINTER(ActionAngleSimulation)]
		self._SimulationStep.restype = None


		self._IntegrateSimulation = self.lib.IntegrateSimulation
		self._IntegrateSimulation.argtypes = [POINTER(ActionAngleSimulation),c_double] 
		self._IntegrateSimulation.restype = c_int
		
		self._mpow2 = self.lib.mpow2
		self._mpow2.argtypes = [c_double,c_int]
		self._mpow2.restype = c_double

		self._mpow = self.lib.mpow
		self._mpow.argtypes = [c_double,c_int]
		self._mpow.restype = c_double
		
		self._secularF2 = self.lib.secularF2
		self._secularF2.argtypes = [c_double]
		self._secularF2.restype = c_double

		self._secularF10 = self.lib.secularF10
		self._secularF10.argtypes = [c_double]
		self._secularF10.restype = c_double

	def secular_coefficients(self,alpha):
		assert alpha < 1.
		try:
			return self._secularF2(alpha),self._secularF10(alpha)
		except:
			print "Failed on input", alpha

		
	def MEGNO_Integration(self,tfin,dt,period,ecc,mu1,mu2,Omega2):
		try:
			return self._MEGNO_Integration(tfin,dt,period,ecc,mu1,mu2,Omega2)
		except:
			print "FAILED ON INPUT: ",tfin,dt,period,ecc,mu1,mu2,Omega2
			return -1.



	def Setup_Integration_Analytic(self,m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,Include0thIn,Include0thOut,resonances1,resonances2,dt):
		Nres1 = resonances1.shape[0];
		Nres2 = resonances2.shape[0];

		arrIn=resonances1.astype(c_int).reshape(-1)
		arrOut=resonances2.astype(c_int).reshape(-1)

		if Include0thIn:
			zflagIn = 1
		else:
			zflagIn = 0
		if Include0thOut:
			zflagOut = 1
		else:
			zflagOut= 0


		sim = ActionAngleSimulation()
		
		varpi2=w2-w1
		varpiTp=wtp-w1

		L0,l0,X0,Y0 = 2.0,lambdatp,np.sqrt(2.)*etp*np.cos(varpiTp),np.sqrt(2.)*etp*np.sin(-1*varpiTp)

		self._InitializeActionAngleSimulation(pointer(sim),Nres1,zflagIn,arrIn,Nres2,zflagOut,arrOut,m1,m2,n1,n2,e1,e2,lambda2,varpi2,L0,l0,X0,Y0,dt)

		return sim


			
	def MEGNO_Integration_Analytic_Full(self,m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,Include0thIn,Include0thOut,resonances1,resonances2,dt,tFin):
		Nres1 = resonances1.shape[0];
		Nres2 = resonances2.shape[0];

		arrIn=resonances1.astype(c_int).reshape(-1)
		arrOut=resonances2.astype(c_int).reshape(-1)
		if Include0thIn:
			zflagIn = 1
		else:
			zflagIn = 0
		if Include0thOut:
			zflagOut = 1
		else:
			zflagOut= 0
		
		sim = ActionAngleSimulation()
		
		varpi2=w2-w1
		varpiTp=wtp-w1
		
		L0,l0,X0,Y0 = 2.0,lambdatp,np.sqrt(2.)*etp*np.cos(varpiTp),np.sqrt(2.)*etp*np.sin(-1*varpiTp)
		self._InitializeActionAngleSimulation(pointer(sim),Nres1,zflagIn,arrIn,Nres2,zflagOut,arrOut,m1,m2,n1,n2,e1,e2,lambda2,varpi2,L0,l0,X0,Y0,dt)
		try:
			self._IntegrateSimulation(pointer(sim),tFin)
			return sim.megno.megno
		except:
			print "FAILED ON INPUT: ",n1,n2,mu1,mu2,resonances1,resonances2,tFin
			return -1.

	def Integrate_Simulation(self,sim,tStop,Nout,MAX_MEGNO=-1):
		megnos=np.zeros(Nout)
		times = np.linspace(0,tStop,Nout)
		for i,t in enumerate(times):
			try:
				self._IntegrateSimulation(pointer(sim),t)
				megnos[i] =  sim.megno.megno
			except:
				print "FAILED ON INPUT: ",n1,n2,mu1,mu2,resonances1,resonances2,tFin
				return -1.
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

	def Integrate_Simulation_GetOrbit(self,sim,tStop,Nout,MAX_MEGNO=-1):
		megnos=np.zeros(Nout)
		orbit = np.zeros((Nout,5))
		times = np.linspace(0,tStop,Nout)
		for i,t in enumerate(times):
			try:
				self._IntegrateSimulation(pointer(sim),t)
				megnos[i] =  sim.megno.megno
				c,s = np.cos(sim.t),np.sin(sim.t)
				R = np.array([[c,-s],[s,c]])
				x,y  = np.dot(R,np.array([sim.state.X,sim.state.Y]))
				orbit[i] = sim.t,sim.state.L,sim.state.l,x,y
			except:
				print "FAILED ON INPUT: ",n1,n2,mu1,mu2,resonances1,resonances2,tFin
				return -1.
			if MAX_MEGNO > 0 and megnos[i] > MAX_MEGNO:
				break
		if i+1==Nout:
			Nfit=int(np.floor(Nout/2))
			tfit=times[Nfit:]
			megnofit=megnos[Nfit:]
			tLy= 1. / np.linalg.lstsq(np.vstack((  tfit,   np.ones(len(tfit))   )).T,megnofit)[0][0]
		else:
			tLy= 1. / np.linalg.lstsq(np.vstack((times[:i+1],np.ones(i+1))).T,megnos[:i+1])[0][0]
		return orbit[:i+1,],megnos[i],tLy

if __name__=="__main__":
	
	w = libwrapper()
	sim=ActionAngleSimulation()


	m1=2.e-5
	m2=2.e-5

	e1=0.01
	e2=0.01
	w1=np.pi/2.

	w2=-np.pi / 4.
	lambda2= w2

	etp=0.02
	lambdatp=0
	wtp=0



	delta1 = 0.013;
	delta2 = 0.01;
	n1 = 3. / 2. * (1+delta1)
	n2 = 2. / 3. / (1+delta2)

	dt=2.*np.pi / 20.

	res1=[]
	res2=[]
	for i,o in Farey_Sequence_N(7):
		j1= o*(3) + i
		res1= res1 + full_resonance_list(j1,o)
		res2= res2 + full_resonance_list(j1,o)
	res1 = np.array(res1,dtype=int)
	res2 = np.array(res2,dtype=int)

	
	sim = w.Setup_Integration_Analytic(m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,True,True,res1,res2,dt)
#	w._SimulationStep(pointer(sim))
	import time
	start_time = time.time()
	data=w.Integrate_Simulation_GetOrbit(sim,2*np.pi*3e4,100,200)
	finish_time = time.time()
	print "<Y>=",data[1],"tLy=",data[2], "time: --- %s seconds ---" %(finish_time - start_time)

# 	orbit = data[0]
# 	plt.plot(orbit[:,0],orbit[:,3])
# 	plt.plot(orbit[:,0],orbit[:,4])
# 	plt.figure()
# 	plt.plot(orbit[:,0],orbit[:,1])
# 	plt.show()
