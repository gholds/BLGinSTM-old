from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments
from Materials import Graphene
import numpy as np
import matplotlib.pyplot as plt

d1, d2 = 1, 305
e1, e2 = 1, 3.9
T = 0
W = 5


BLG = Graphene.Bilayer()
exp1 = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,W,screening=False)
exp2 = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,W,screening=True)

vplus = 0.1
vminus = 0.01
VT = np.array([-0.05,0.05])
VB = np.array([5,10])

exp1.generate_tunnelcurrent(VT,VB)
exp2.generate_tunnelcurrent(VT,VB)

print(exp1.I)
print(exp2.I)
