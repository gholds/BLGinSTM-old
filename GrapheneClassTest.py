from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments
from StatisticalDistributions import Temperature
from UniversalConstants import *

import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

BLG = Graphene.Bilayer()

n =10**16

# charge = lambda x: BLG.nplusT0(x,vminus) - n
# vplus = optimize.newton(charge,0.1)
# print(vplus)
# print(BLG.screened_newton(vplus,vminus))

VMext = np.linspace(-0.1,0.1,num=100)

vm1 = []
vm2 = []

# for vmext in VMext:
# 	# Find value of vplus that yields this carrier density
# 	vm = BLG.screened_vminus(n,vmext)
# 	vm1.append(vm)
# 	# vm1.append(Temperature.FermiDirac(KE-q*abs(vplus),0))

VMext = np.array([0.1])
vm1 = BLG.screened_vminus(n,VMext)
#vm2 = BLG.screened_vminus(5*n,VMext)

print(vm1[0])

# plt.plot(VMext,vm1)
# plt.plot(VMext,vm2)
# plt.text(-0.07,0.08,'n={:.2E}'.format(n))

# plt.show()