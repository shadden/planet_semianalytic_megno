from ctypes import *
import numpy as np

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
	w = libwrapper()
	n1range = np.linspace(1.5*(1-0.007),1.5*(1+0.007),10)
	n2range = np.linspace(1./1.5/(1+0.007),1./1.5/(1-0.007),10)

	tfin = 2*np.pi*1e+4
	dt = 2*np.pi/30.
	m1=1e-5
	m2=1e-5
	
	def f(pars):
		n1,n2 = pars
		period = n1
		omega2 = n2 / n1
		return w.MEGNO_Integration(tfin,dt,period,0.,m1,m2,omega2)
	
	from rebound.interruptible_pool import InterruptiblePool
	pool = InterruptiblePool()
	parameters = [ (n1,n2) for n1 in n1range for n2 in n2range ]
	results = pool.map(f , parameters)

	resultsArr = np.hstack((np.array(parameters),np.array(results).reshape(-1,1) ))
	resultsArr.tofile("results.dat")


