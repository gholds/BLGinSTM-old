import numpy as np

# Various directories
cd_folder = 'data/CarrierDensity/'

def get_vplus(T):
    vplus = np.load('data/CarrierDensity/Temperature_{:.1E}K/vplus.npy'.format(T))
    return np.concatenate((-vplus[:0:-1],vplus))

def get_vminus(T):
    vminus = np.load('data/CarrierDensity/Temperature_{:.1E}K/vminus.npy'.format(0))
    return np.concatenate((-vminus[:0:-1],vminus))

def get_nplus(T):
    nplus_surface = np.load('data/CarrierDensity/Temperature_{:.1E}K/nplus_surface.npy'.format(0))
    nplus_surface = np.concatenate((nplus_surface[:0:-1,:],nplus_surface))
    nplus_surface = np.concatenate((-nplus_surface[:,:0:-1],nplus_surface),axis = 1)
    return nplus_surface

def get_nminus(T):
    nminus_surface = np.load('data/CarrierDensity/Temperature_{:.1E}K/nminus_surface.npy'.format(0))
    nminus_surface = np.concatenate((-nminus_surface[:0:-1,:],nminus_surface))
    nminus_surface = np.concatenate((nminus_surface[:,:0:-1],nminus_surface),axis = 1)
    return nminus_surface