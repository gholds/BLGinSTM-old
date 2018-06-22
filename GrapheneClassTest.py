from Materials import Graphene
from TunnelingExperiment import TunnelingExperiments

import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

BLG = Graphene.Bilayer()

n = 1*10**16

# charge = lambda x: BLG.nplusT0(x,vminus) - n
# vplus = optimize.newton(charge,0.1)
# print(vplus)
# print(BLG.screened_newton(vplus,vminus))

VMext = np.linspace(-0.1,0.1,num=100)

vm1 = []
vm2 = []

for vmext in VMext:
	# Find value of vplus that yields this carrier density
	charge = lambda x: BLG.nplusT0(x,vmext) - n
	vplus = optimize.newton(charge,0.1)

	print(vplus,vmext)
	vm1.append(BLG.screened_vminus(vplus,vmext))
	#vm2.append(BLG.screened_newton(vplus,vmext))

vm1 = np.array(vm1)

plt.plot(VMext,vm1,label='Abergeletc')
plt.text(-0.07,0.08,'n={:.2E}'.format(n))
#plt.plot(VMext,vm2,label='Newton')

plt.show()


# from scipy import optimize

# def ftest(x,y):
# 	return (x-1) + (y-2)**2

# arg = (2,)
# print(optimize.newton(ftest,0.5,args=arg))

# d1, d2 = 1, 305
# e1, e2 = 1, 3.9
# T = 0
# Wtip = 5

# stm = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)


# Ks = np.linspace(0,10**9,num=100)
# vminus = 0.05


# pdiff_g		= [] 
# pdiff_LE	= []
# pdiff_none	= []

# for k in Ks:
# 	pdiff_g.append(BLG.Pdiff(k,vminus,approx='g3=0'))
# 	pdiff_LE.append(BLG.Pdiff(k,vminus,approx='LowEnergy'))
# 	pdiff_none.append(BLG.Pdiff(k,vminus,approx='None'))

# pdiff_g = np.array(pdiff_g)
# pdiff_LE = np.array(pdiff_LE)
# pdiff_none = np.array(pdiff_none)

# plt.plot(Ks,pdiff_g,label='Current')
# plt.plot(Ks,pdiff_LE,label='Young and Levitov')
# plt.plot(Ks,pdiff_none,label='None')
# plt.title('Probability Differences')
# plt.xlabel('k (1/m)')
# plt.ylabel('Prob Difference')
# plt.legend(title='Approximation')
# plt.show()

# d1 	= 1
# d2 	= 305
# e1 	= 1
# e2 	= 3.9
# T 	= 0
# Wtip= 5



# stm = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)


# VTrange = [-0.2,0.2]
# num_vts_100 = 1

# VBrange = [-30,30]
# num_vbs_100 = 1

# results = stm.generate_tunnelcurrent(VTrange,num_vts_100,VBrange,num_vbs_100)


# plt.imshow(np.gradient(results,axis=1))
# plt.show()
