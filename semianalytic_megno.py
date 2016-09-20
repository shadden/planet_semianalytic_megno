from ctypes import *
import numpy as np

import os
who =os.popen("whoami") 
if who.readline().strip() =='samuelhadden':
	SRCDIR = "/Users/samuelhadden/02_ThreePlanet_Stability/03_Semianalytic_Code"
else:
	SRCDIR = "/Users/samuelhadden/02_ThreePlanet_Stability/03_Semianalytic_Code"
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
		
		
