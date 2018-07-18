from TunnelingExperiment import TunnelingExperiments
import numpy as np
d1, d2 = 0.5, 305
e1, e2 = 1, 3.9
T = 0
Wtip = 5

experiment = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)

VT = np.linspace(-0.2,0.2, num=20)

VB = np.linspace(-40,40,num=20)

experiment.generate_tunnelcurrent(VT,VB)

experiment.plot_dIdV()