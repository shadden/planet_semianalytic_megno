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
				("varpi2",c_double)]
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

	res1= np.array([[3,1,0]]) #np.array([[3,1,0],[3,1,1],[6,2,0],[6,2,1],[6,2,2]])
	res2=np.array([[3,1,0],[3,1,1],[6,2,0],[6,2,1],[6,2,2]])
	

	
	m1=0.e-5
	m2=3.e-6

	e1=0.0
	e2=0.04
	w1=w2=0.
	lambda2= 0.

	etp=0.02
	lambdatp=3.5
	wtp=2.


	
	delta1 = -0.001;
	delta2 = 0.005;
	n1 = 3. / 2. * (1+delta1)
	n2 = 2. / 3. / (1+delta2)

	

	dt=2.*np.pi / 50.
	tFin = 2*np.pi * 1500. *2   
	Nsteps = int( tFin / dt )
	Ndump = 3*20 * 10 
	
	import rebound
	def SimulationSetup():
		sim = rebound.Simulation()
		sim.units= ('yr','AU','Msun')
		sim.integrator = "whfast"
		sim.integrator_whfast_safe_mode = 0
		sim.exit_max_distance = 3.
		sim.dt = dt / 2. / np.pi
		sim.add(m=1.0);
		sim.add(m=m1,id=1,a=n1**(-2./3.),e=e1,l=0);
		sim.add(m=m2,id=2,a=n2**(-2./3.),e=e2,l=lambda2);
		sim.add(m=0.,id=3,a=1.,e=etp,pomega=wtp,l=lambdatp);
		sim.move_to_com()
		return sim



	sim = w.Setup_Integration_Analytic(m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,True,True,res1,res2,dt,tFin)
	simNbody = SimulationSetup()
		
	
	Npts=int(np.floor(Nsteps/Ndump)) 
	data = np.zeros((4,Npts))
	NBdata = np.zeros((3,Npts))
	j=0	
	for i in range(Nsteps):
		if i%Ndump==0 and j<Npts:

			data[0,j] = i*dt
			data[1,j] = sim.state.L
			data[2,j] = np.sqrt(sim.state.X*sim.state.X + sim.state.Y*sim.state.Y) / np.sqrt(2)
			data[3,j] = sim.megno.megno
			orbs=simNbody.calculate_orbits()
			NBdata[0,j] = i*dt
			NBdata[1,j] =orbs[-1].a
			NBdata[2,j] =orbs[-1].e
			simNbody.integrate(i*dt/2./np.pi,exact_finish_time=0)
			j=j+1
			


			
		w._SimulationStep(pointer(sim))

	fig,ax = plt.subplots(3,1,sharex=True)	
	a0 = NBdata[1,0]
	sma = a0 * 0.25 * (data[1])**2
	ax[0].plot(NBdata[0],NBdata[1],'k-')
	ax[0].plot(data[0],sma,'r-')
# 
	ax[1].plot(NBdata[0],NBdata[2],'k-')
	ax[1].plot(data[0],data[2],'r-')
# 
	ax[2].plot(data[0],data[3],'r-')
	ax[2].set_ylim((1.5,3.))

	plt.show()
	

	
if False: #__name__=="__main__":
	
	w = libwrapper()
	sim=ActionAngleSimulation()

	Ngrid=15
	n10 =1.5
	n20 =2./3.

	res1=np.array([[3,1,0],[3,1,1]])
	res2=np.array([[3,1,0],[3,1,1]])
	
	dw = 1.74532925199
	
	m1=1.e-5
	m2=1.e-5
	e0=e1=e2=0.04

	lambda2= 2 * dw
	varpi2= 2 * dw
	

	dt=2.*np.pi / 20.
	
	pars=[]
	for delta2 in np.linspace(+0.02,-0.02,Ngrid):
		for delta1 in np.linspace(-0.02,0.02,Ngrid):
			n1 = n10 * (1+delta1)
			n2 = n20 / (1+delta2)
			pars.append((n1,n2))

# 	def MEGNO_Integration_Analytic_Full(self,m1,m2,n1,n2,e1,e2,etp,w1,w2,wtp,lambda2,lambdatp,Include0thIn,Include0thOut,resonances1,resonances2,dt,tFin):
	def f(x):
		n1,n2=x
		meg=w.MEGNO_Integration_Analytic_Full(m1,m2,n1,n2,e1,e2,e0,0.,varpi2,dw,lambda2,dw,True,True,res1,res2,dt,2*np.pi*3e3)
		return meg
	
	from rebound.interruptible_pool import InterruptiblePool
	pool = InterruptiblePool()
	results=pool.map(f,pars)
	results2d = np.array(results).reshape(Ngrid,Ngrid)

	fig = plt.figure(figsize=(7,5))
	ax = plt.subplot(111)
	delta=0.01
	extent = [ n10 * (1-delta), n10 * (1+delta), n20 / (1+delta), n20 / (1-delta)]
	ax.set_xlim(extent[0],extent[1])
	ax.set_xlabel("inner freq $n1$")
	ax.set_ylim(extent[2],extent[3])
	ax.set_ylabel("outer freq $n2$")
	im = ax.imshow(results2d, interpolation="none", vmin=1.9, vmax=4, cmap="RdYlGn_r", origin="lower", aspect='auto', extent=extent)
	cb = plt.colorbar(im, ax=ax)
	cb.set_label("MEGNO $\\langle Y \\rangle$")
	plt.show()
	
