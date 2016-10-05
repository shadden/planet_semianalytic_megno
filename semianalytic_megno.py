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
	_fields_ = [("Nres",c_int),
				("MAxOrder",c_int),
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
		self._InitializeActionAngleSimulation.argtypes = [POINTER(ActionAngleSimulation),c_int,p1dInt,c_int,p1dInt] + [c_double for i in range(12)]
		self._InitializeActionAngleSimulation.restype = None

		self._SimulationStep = self.lib.SimulationStep
		self._SimulationStep.argtypes = [POINTER(ActionAngleSimulation)]
		self._SimulationStep.restype = None


		self._IntegrateSimulation = self.lib.IntegrateSimulation
		self._IntegrateSimulation.argtypes = [POINTER(ActionAngleSimulation),c_double] 
		self._IntegrateSimulation.restype = c_int
		
	def MEGNO_Integration(self,tfin,dt,period,ecc,mu1,mu2,Omega2):
		try:
			return self._MEGNO_Integration(tfin,dt,period,ecc,mu1,mu2,Omega2)
		except:
			print "FAILED ON INPUT: ",tfin,dt,period,ecc,mu1,mu2,Omega2
			return -1.

	def MEGNO_Integration_Analytic(self,n1,n2,m1,m2,resonances1,resonances2,dt,tFin):
		Nres1 = resonances1.shape[0];
		Nres2 = resonances2.shape[0];

		arrIn=resonances1.astype(c_int).reshape(-1)
		arrOut=resonances2.astype(c_int).reshape(-1)
		
		sim = ActionAngleSimulation()
		
		L0,l0,X0,Y0 = 2.0,0.,0.0,0.0
		e1=0
		e2=0
		varpi2=0
		
		self._InitializeActionAngleSimulation(pointer(sim),Nres1,arrIn,Nres2,arrOut,m1,m2,n1,n2,e1,e2,varpi2,L0,l0,X0,Y0,dt)
		try:
			self._IntegrateSimulation(pointer(sim),tFin)
			return sim.megno.megno
		except:
			print "FAILED ON INPUT: ",n1,n2,mu1,mu2,resonances1,resonances2,tFin
			return -1.
		
if __name__=="__main__":
	
	w = libwrapper()
	sim=ActionAngleSimulation()
	res1=np.vstack((np.arange(4,7),np.ones(3),np.zeros(3))).T
	res2=np.vstack((np.arange(4,7),np.ones(3),np.ones(3))).T
	m1=1e-5
	m2=3.e-5
	e1=e2=0
	varpi2=0
	L0=2
	X0=Y0=l0=0
	n1=5./4.*(1+0.005)
	n2=4./5.*(1+0.009)
	dt=2.*np.pi / 10.
	
	results=[]
	for delta in np.linspace(-0.005,0.005,20):
		n1=5./4.*(1+delta)
		meg=w.MEGNO_Integration_Analytic(n1,n2,m1,m2,res1,res2,dt,2*np.pi*1e4)
		print n1,meg
		results.append((n1,meg))

	results=np.array(results)
	plt.plot(results[:,0],results[:,1])
	plt.show()
	
if False: #__name__=="__main__":

	parser = ArgumentParser(description='Run a grid of simulations and compute the MEGNO')
	parser.add_argument('-N','--Ngrid',metavar='N',type=int,default=10,help='Number of grid points in each dimension')
	parser.add_argument('-C','--checkpoints',metavar='N',type=int,default=0,help='Number of checkpoints')
	parser.add_argument('--restart', default=False, action='store_true', help='continue a previously-existing run')
	parser.add_argument('-T','--Time',metavar='t',type=float,default=5.0e4,help='Lengh of simulations')
	parser.add_argument('-F','--infile',metavar='<FILE>',default='input.txt',help='Text file with period ranges')

	args = parser.parse_args()
	Ngrid = args.Ngrid
	simLength = args.Time * 2 * np.pi
	infile = args.infile
	checkpoints = args.checkpoints + 1
	restart = args.restart

	with open(infile,'r') as fi:
		planetMass1,planetMass2 = map(float, fi.readline().split() )
		n1min,n1max = map(float, fi.readline().split() )
		n2min,n2max = map(float, fi.readline().split() )
		e1,w1  = map(float, fi.readline().split() )
		e2,w2  = map(float, fi.readline().split() )
		etp,wtp  = map(float, fi.readline().split() )
		theta1,theta2,thetatp=map(float,fi.readline().split())
		Mtp  = map(float, fi.readline().split() )[0]

	print "Sim. length:", simLength
	print "Planet 1: m=%.1e , e=%.3f , w= %.3f deg. , theta= %.3f deg."%(planetMass1,e1,w1,theta1)
	print "Planet 2: m=%.1e , e=%.3f , w= %.3f deg. , theta= %.3f deg."%(planetMass2,e2,w2,theta2)
	print "Test particle: m=%.1e , e=%.3f , w= %.3f deg. , theta= %.3f deg."%(Mtp,etp,wtp,thetatp)

	w1 = np.pi * w1 / 180.
	w2 = np.pi * w2 / 180.
	wtp = np.pi * wtp / 180.
	theta1 = np.pi * theta1 / 180.
	theta2 = np.pi * theta2 / 180.
	thetatp = np.pi * thetatp / 180.


	rHill = ( np.max([ planetMass1,planetMass2 ]) /3.)**(1./3.)

	w = libwrapper()

# 	def f(pars):
# 		n1,n2 = pars
# 		n1 = n1 / 2. / np.pi
# 		n2 = n2 / 2. / np.pi
# 		period = n1
# 		omega2 = n2 / n1
# 		return w.MEGNO_Integration(simLength, 2*np.pi/30. ,period,0.,planetMass1,planetMass2,omega2)
	
	def f(pars):
		n1,n2 = pars
		n1 = n1 / 2. / np.pi
		n2 = n2 / 2. / np.pi
		return w.MEGNO_Integration_Analytic(n1,n2,planetMass1,planetMass2,np.array([7,8],dtype=c_int),np.array([8,9],dtype=c_int),simLength)


	par_d1 =  np.linspace(n1min,n1max,Ngrid+1)[:-1] 
	par_d2 =  np.linspace(n2min,n2max,Ngrid+1)[:-1]

	parameters = []
	for d2 in par_d2:
		for d1 in par_d1:
			parameters.append((d1,d2))

	parameters = np.array(parameters)
	from rebound.interruptible_pool import InterruptiblePool	
	pool = InterruptiblePool()
	if restart:
		import re
		import glob
		checkpointfiles = glob.glob("checkpoint_*.dat")
		exprn = re.compile("checkpoint_(\d+)")
		checkPointNumbers = [ int(re.search(exprn,x).group(1)) for x in checkpointfiles ]
		cMax = np.max(checkPointNumbers)
		results = np.fromfile("checkpoint_%d.dat"%cMax).tolist()
		print "Restarting from checkpoint file checkpoint_%d.dat..."%cMax
		sys.stdout.flush()
	else:
		results = []
		cMax = -1

	for c in range(checkpoints)[cMax+1:]:
		index_low = c*len(parameters)/checkpoints
		index_high = (c+1)*len(parameters)/checkpoints
		results = results + pool.map(f , parameters[index_low:index_high])
		np.array(results).tofile("checkpoint_%d.dat"%c)
	
	resultsArr = np.hstack((np.array(parameters),np.array(results).reshape(-1,1) ))
	np.save("MEGNO_RESULTS.npy",resultsArr)
	resultsArr.tofile("MEGNO_RESULTS.dat")
