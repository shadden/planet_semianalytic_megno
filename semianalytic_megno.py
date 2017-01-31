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
	_fields_ = [("state", ActionAnglePhaseState),
				("parameters",SimulationParameters),
				("megno",MEGNO_Auxilary_Variables),
				("rIn",ResonanceData),
				("rOut",ResonanceData),
				("dt",c_double),
				("t",c_double)]


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



	def Setup_Integration_Analytic(self,m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,Include0thIn,Include0thOut,resonances1,resonances2,dt,tFin):
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


if __name__=="__main__":
	
	w = libwrapper()
	sim=ActionAngleSimulation()

# 	res1= np.array([[3,1,0],[3,1,1]]) #np.array([[3,1,0],[3,1,1],[6,2,0],[6,2,1],[6,2,2]])
# 	res2=np.array([[3,1,0],[3,1,1],[6,2,0],[6,2,1],[6,2,2]]) #
	res1 =np.array([])	
	res2 =np.array([])	

	
	m1=0.e-6
	m2=1.e-5

	e1=0.0
	e2=0.05
	w1=0.
	w2=-np.pi / 4.
	lambda2= w2

	etp=0.02
	lambdatp=0
	wtp=0



	delta1 = -0.001;
	delta2 = 0.01;
	n1 = 3. / 2. * (1+delta1)
	n2 =  1/2.15  #2. / 3. / (1+delta2)

	dt=2.*np.pi / 20.
	tFin = 2*np.pi * 1500. * 150   

	Nsteps = int( tFin / dt )
	Npts = 2**10
	Ndump = int(np.floor(Nsteps / Npts))


	sim = w.Setup_Integration_Analytic(m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,True,True,res1,res2,dt,tFin)
	w._SimulationStep(pointer(sim))
	print sim.state.X,sim.state.Y
	print sim.parameters.lambda2 ,  sim.parameters.varpi2
# 	import rebound
# 	def SimulationSetup():
# 		sim = rebound.Simulation()
# 		sim.units= ('yr','AU','Msun')
# 		sim.integrator = "whfast"
# 		sim.integrator_whfast_safe_mode = 0
# 		sim.exit_max_distance = 3.
# 		sim.dt = dt / 2. / np.pi
# 		sim.add(m=1.0);
# 		sim.add(m=m1,id=1,a=n1**(-2./3.),e=e1,l=0);
# 		sim.add(m=m2,id=2,a=n2**(-2./3.),e=e2,pomega=w2,l=lambda2);
# 		sim.add(m=0.,id=3,a=1.,e=etp,pomega=wtp,l=lambdatp);
# 		sim.move_to_com()
# 		return sim
# 
# 
# 
# 	sim = w.Setup_Integration_Analytic(m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,True,True,res1,res2,dt,tFin)
# 	simNbody = SimulationSetup()
# 		
# 	data = np.zeros((4,Npts))
# 	NBdata = np.zeros((3,4,Npts))
# 	j=0	
# 
# 	for i in range(Nsteps):
# 		if i%Ndump==0 and j<Npts:
# 
# 			data[0,j] = i*dt
# 			data[1,j] = sim.state.L
# 			data[2,j] = np.sqrt(sim.state.X*sim.state.X + sim.state.Y*sim.state.Y) / np.sqrt(2)
# 			data[3,j] =  i * n1 * dt - np.arctan2(sim.state.Y,sim.state.X) 
# 	
# 			orbs=simNbody.calculate_orbits()
# 			for k in range(3):
# 				NBdata[k,0,j] = i*dt
# 				NBdata[k,1,j] =orbs[k].a
# 				NBdata[k,2,j] =orbs[k].e
# 				NBdata[k,3,j] =orbs[k].pomega
# 	
# 			simNbody.integrate(i*dt/2./np.pi,exact_finish_time=0)
# 			j=j+1
# 		w._SimulationStep(pointer(sim))
# 
# 	fig,ax = plt.subplots(2,1,figsize=(10,5))
# 	ax[0].plot(NBdata[2,0],NBdata[2,2]*np.cos(NBdata[2,3]),'k-')
# 	ax[0].plot(NBdata[2,0],NBdata[2,2]*np.sin(NBdata[2,3]),'r-')	
# 	ax[0].plot(data[0],data[2]*np.cos(data[3]),'b-')
# 	ax[0].plot(data[0],data[2]*np.sin(data[3]),'g-')
# 	for i in range(2):
# 		if i==0:
# 			style = '-'
# 		else:
# 			style = "--"
# 		ax[1].plot(NBdata[i,0],NBdata[i,2]*np.cos(NBdata[i,3]),c='k',ls=style)
# 		ax[1].plot(NBdata[i,0],NBdata[i,2]*np.sin(NBdata[i,3]),c='r',ls=style)	
# 	plt.show()
	
