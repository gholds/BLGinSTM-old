from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


BLG = Graphene.Bilayer()



d1 	= 1
d2 	= 305
e1 	= 1
e2 	= 3.9
T 	= 0
Wtip= 5



stm = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)


VTrange = [-0.2,0.2]
num_vts_100 = 1

VBrange = [-30,30]
num_vbs_100 = 1

results = stm.generate_tunnelcurrent(VTrange,num_vts_100,VBrange,num_vbs_100)


plt.imshow(np.gradient(results))
plt.show()
