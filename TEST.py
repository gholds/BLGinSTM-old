from Materials import Graphene
import numpy as np
import matplotlib.pyplot as plt

BLG = Graphene.Bilayer()

n = 10**16

vminus = np.linspace(-0.1,0.1,num=100)

vminus_screened = BLG.screened_vminus(n,vminus)

plt.plot(vminus, vminus_screened)

plt.show()