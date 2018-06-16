from Universal_Constants import *
from BandStructure import *

u = 0*eVtoJ
print('Abs Min: ', round(emin(u)*JtoeV,3))
print('Loc Max: ', u/2 * JtoeV)

vplus = 0.03

print('Fermi Level: ',vplus)

print(carrierDensity(q*vplus,u))
print(nplusT0(q*vplus,u))