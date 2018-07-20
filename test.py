from TunnelingExperiment import STM
import numpy as np
d1, d2 = 0.5, 305
e1, e2 = 1, 3.9
T = 0
Wtip = 5

experiment = STM.BLGinSTM(d1,d2,e1,e2,T,Wtip)

VT = np.linspace(-0.5,0.5, num=100)
VB = np.linspace(-40,40,num=10)

experiment.plot_dIdV_waterfall(VT,VB)
#experiment.generate_tunnelcurrent(VT,VB)
#experiment.plot_dIdV()