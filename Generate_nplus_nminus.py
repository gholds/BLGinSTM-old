# This file is a script to generate and save the surfaces defined by nplus and nminus
import numpy as np
from BLG.BandStructure import nplus, nminus, nplusT0
from BLG.BLG_Constants import q

def generate_nplus_nminus(T):
	''' Generates and stores the nplus and nminus surfaces at
	Temperature T'''
	pass

# Directory for all Carrier Density Surfaces
cd_folder = 'data/CarrierDensity/'

# Due to symmetry and antisymmetry, only need the first quadrant
vplus = np.linspace(0,0.4,num=800)
vminus = np.linspace(0,0.25,num=500)[:,np.newaxis]

# Set temperature
T=0

# Create a directory if it does not exist
import os

# Create folder to save carrier density surfaces
temp_folder = cd_folder + 'Temperature_{:.1E}K/'.format(T)

# Check that the directory exists
if not os.path.exists(temp_folder):
	os.makedirs(temp_folder)



# Choose the size of the batches we will generate
d = 10

# Check that it is compatible with the lengths of vplus and vminus
if len(vplus) % d != 0 or len(vminus) % d != 0:
	print('Batch size (d) incompatible with voltage arrays')
	print('d= {} does not evenly divide len(vplus)= {} or len(vminus) = {}'.format(d,len(vplus),len(vminus)))
	import sys
	sys.exit()


# T=0 is special because we can use nplusT0
if T == 0:
	# nplusT0 is capable of broadcasting and is fast
	nplus_surface = nplusT0(vplus,vminus)

	# nminus is also capable of broadcasting, but is
	# not fast. We use a for-loop
	nminus_surface = np.empty(np.shape(vminus*vplus))

	for i in range(int(len(vplus)/d)):
		for j in range(int(len(vminus)/d)):
			nminus_surface[d*j:d*j+d,d*i:d*i+d]=nminus(vplus[d*i:d*i+d],vminus[d*j:d*j+d,:],0)

# Otherwise, temperature nonzero, and all integrals numerical
else:

	nplus_surface = np.empty(np.shape(vminus*vplus))
	nminus_surface = np.empty(np.shape(vminus*vplus))

	for i in range(int(len(vplus)/d)):
		for j in range(int(len(vminus)/d)):
			nplus_surface[d*j:d*j+d,d*i:d*i+d]=nplus(vplus[d*i:d*i+d],vminus[d*j:d*j+d,:],T)
			nminus_surface[d*j:d*j+d,d*i:d*i+d]=nminus(vplus[d*i:d*i+d],vminus[d*j:d*j+d,:],T)



# Save the surfaces
np.save(temp_folder+'nplus_surface.npy',nplus_surface)
np.save(temp_folder+'nminus_surface.npy',nminus_surface)

# Save the voltages
np.save(temp_folder+'vplus.npy', vplus)
np.save(temp_folder+'vminus.npy', vminus)

