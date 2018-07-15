from TunnelingExperiment import TunnelingExperiments
import numpy as np

d1, d2 = 1, 305
e1, e2 = 1, 3.9
T = 0
W = 5

experiment = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,W)

VT = np.linspace(-1,1, num=10)
VB = np.linspace(-45,45,num=10)

vs = experiment.generate_vplus_vminus([VT[0],VT[-1]],0.05,[VB[0],VB[-1]],0.05,method='DasSarma')
print(vs[0])