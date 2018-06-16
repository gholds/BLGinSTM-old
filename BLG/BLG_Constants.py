import numpy as np
from BLG.Universal_Constants import *

a=1.42 * (10**(-10)) # m, Interatom spacing
Ac = 3*np.sqrt(3)*(a**2) / 2 # Area of unit cell of graphene
g0=2.8*eVtoJ # J, Interatom hopping potential

vF=3*a*g0/(2*hbar) # m/s, Fermi velocity

g1 = 0.358*eVtoJ # J, A1-A2 hopping
a0 = 3*10**(-10) # m, distance between graphene layers

randomnumber = 444