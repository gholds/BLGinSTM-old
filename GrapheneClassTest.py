from Materials import Graphene

import numpy as np
from BLG.Universal_Constants import eVtoJ

blg = Graphene.Bilayer()

vplus = np.linspace(0,0.4,num=100)
vminus = 0.005 

nplus_g30 = blg.nplus(vplus,vminus,0)
nplus_lowe = blg.nplus(vplus,vminus,0,approx='LowEnergy')

nminus_g30 = blg.nminus(vplus,vminus,0)
nminus_lowe = blg.nminus(vplus,vminus,0,approx='LowEnergy')

import matplotlib.pyplot as plt

plt.plot(vplus,nplus_g30, label = 'nplus g3=0')
plt.plot(vplus,nplus_lowe, label = 'nplus, LowEnergy')

plt.plot(vplus,nminus_g30, label = 'nminus g3=0')
plt.plot(vplus,nminus_lowe, label = 'nminus Low E')

plt.legend()
plt.show()