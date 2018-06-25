from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments
from StatisticalDistributions import Temperature
from UniversalConstants import *

import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize


d1, d2 = 1, 305 # nanometers
e1, e2 = 1, 3.9 # relative permittivity
T = 0			# K
Wtip = 5 		# eV

STM = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)

points = 30
points = points/100

I =STM.generate_tunnelcurrent([-0.2,0.2],points,[-30,30],points,method='DasSarma')

plt.imshow(np.gradient(I,axis=0))

plt.show()