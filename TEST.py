from Materials import Graphene
import numpy as np
import matplotlib.pyplot as plt

BLG = Graphene.Bilayer()

n1 = 10**16
n2 = 5 * n1

vminus = np.linspace(-0.1,0.1,num=100)

vms1 = BLG.screened_vminus(n1,vminus)
vms2 = BLG.screened_vminus(n2,vminus)

plt.plot(vminus, vms1)
plt.plot(vminus, vms2)

plt.show()