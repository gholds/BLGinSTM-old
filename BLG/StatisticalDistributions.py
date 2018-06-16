import numpy as np

kB = 1.38064852 * 10**-23 # Boltmann constant J/K


def FermiDirac(e,T):
    '''Returns the Fermi-Dirac distribution.'''
    # Using 1 / (np.exp()+1) introduces overflow errors. This keeps the floats of managable size.
    return np.exp(-np.logaddexp(e/(kB*(T+0.0001)),0))
