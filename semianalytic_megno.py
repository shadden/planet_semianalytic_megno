from ctypes import *
import numpy as np
from argparse import ArgumentParser
import os
who =os.popen("whoami") 
if who.readline().strip() =='samuelhadden':
	SRCDIR = "/Users/samuelhadden/02_ThreePlanet_Stability/03_Semianalytic_Code"
else:
	SRCDIR = "/projects/p20783/sjh890/02_Chaos_Project/SemiAnalyticCodes/planet_semianalytic_megno"
who.close()


class libwrapper(object):

	def __init__(self):
		
		self.lib = CDLL("%s/libsemianalyticMEGNO.so"%SRCDIR)
		self._MEGNO_Integration = self.lib.MEGNO_Integration
		self._MEGNO_Integration.argtypes = [c_double for i in range(7)]
		self._MEGNO_Integration.restype = c_double
		
	def MEGNO_Integration(self,tfin,dt,period,ecc,mu1,mu2,Omega2):
		try:
			return self._MEGNO_Integration(tfin,dt,period,ecc,mu1,mu2,Omega2)
		except:
			print "FAILED ON INPUT: ",tfin,dt,period,ecc,mu1,mu2,Omega2
			return -1.
		
		
if __name__=="__main__":

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
	def f(pars):
		n1,n2 = pars
		n1 = n1 / 2. / np.pi
		n2 = n2 / 2. / np.pi
		period = n1
		omega2 = n2 / n1
		return w.MEGNO_Integration(simLength, 2*np.pi/30. ,period,0.,planetMass1,planetMass2,omega2)


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
