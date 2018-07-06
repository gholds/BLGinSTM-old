from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments
from StatisticalDistributions import Temperature
from UniversalConstants import *

import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize,integrate


d1, d2 = 1, 305 # nanometers
e1, e2 = 1, 3.9 # relative permittivity
T = 0			# K
Wtip = 5 		# eV

BLG = Graphene.Bilayer()
STM = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)

VT = np.linspace(-0.2,0.2,num=100)

STM.generate_tunnelcurrent([-0.2,0.2],1,
							[30,30],.01,
							method='DasSarma')

plt.plot(VT,np.gradient(STM.I))
plt.show()