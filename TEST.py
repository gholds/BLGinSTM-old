from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments
import numpy as np
import matplotlib.pyplot as plt

d1, d2 = 1, 305
e1, e2 = 1, 3.9
T = 0
W = 5

exp = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,W)

vplus = 0.01
VT = 1
VB = 30

print(exp.vminus_n1(vplus,VT,VB))