import numpy as np
from scipy import optimize

from BLG.BLG_Constants import *
from BLG.Universal_Constants import *

from BLG.StatisticalDistributions import *

from scipy import integrate

A = 30*10**-12

#########################
### Bands and extrema ###
#########################
def Dispersion(k,u,band):
    ''' Returns the energy (J) of an electron with wavevector k (rad/m) in band = 1 or 2.
    To get holes, simply result multiply by -1.'''
    p = hbar * vF * k # define a sort of normalized momentum 
    radical=(g1**4)/4 + (u**2 + g1**2)*(p**2)
    return np.sqrt( (g1**2)/2 + (u**2)/4 + p**2 + ((-1)**(band))*np.sqrt(radical) )

def kmin(u, band=1):
    '''Returns positive wavenumber at the minimum of the first band in 1/m'''
    k2 = ( u**2 / (2*hbar*vF)**2 ) * ( (2*g1**2 + u**2) /( g1**2 + u**2 ) )
    return np.sqrt(k2)

def emin(u):
    '''Returns minimum of the first band in Joules'''
    emin2 = (u/2)**2 * (g1**2 / ( g1**2 + u**2 ) )
    return np.sqrt(emin2)

##################################
### Fermi wavevector and level ###
##################################
def kFermi(n,u,pm):
    ''' Returns Fermi vector kF+ for pm=1 and kF- for pm=2 in units rad/m'''
        
    # Define the more complicated factors and terms
    numerator = (pi * hbar**2 * vF**2 * n)**2 + ( g1*u )**2
    denominator = g1**2 + u**2
    pmterm = 2*pi*hbar**2 * vF**2 * abs(n) # abs for fact that electrons and holes symmetric
    
    # Factor proportional to k**2
    propk2 = ( numerator / denominator ) + u**2 + (-1)**(pm-1) * pmterm
    
    # If the fermi level is above u/2, set kF- to zero
    # This says that the region of occupied states is now a disk
    if pm%2==0:
        propk2 = (propk2 >= 0) * propk2
        propk2 = (Dispersion(kFermi(n,u,1),u,1)<u/2) * propk2
    
    return np.sqrt( propk2 ) / (2*hbar*vF)

def eFermi(n,u):
    '''Returns the Fermi level (Joules) given density n and interlayer potential energy difference u
    Positive n returns a positive Fermi level, meaning positive carrier densities are electrons by convention.'''
    
    numerator = (hbar**2 * vF**2 * n *pi)**2 + (g1 * u)**2
    denominator = 4 * (g1**2 + u**2)
    
    return np.sign(n) * np.sqrt( numerator / denominator )

#######################
### Carrier Density ###
#######################
def carrierDensity(eF,u):
    '''Returns the carrier density per unit area.
    Returns a positive value for positive fermi levels and vice versa.'''
    
    # Calculate the radical
    radical = 4*(g1**2+u**2)*eF**2 - g1**2*u**2
    
    # For energies within the gap, radical is negative, so set it to 1 instead
    radical = (radical>=0)*radical + (radical<0)
    
    # Return also multiplying by the truth value for energies being outside the gap
    return (abs(eF)>emin(u))*np.sign(eF)*(1 / (hbar**2 * vF**2 * pi)) * np.sqrt(radical)


def nplusT0(vplus,vminus):

    eF = q*vplus
    u  = -2*q*vminus
    # Calculate the radical
    radical = (g1**2+u**2) * eF**2 - g1**2 * u**2 / 4

    # For energies within the gap, radical is negative, so set it to 0 instead
    radical = (radical>=0)*radical

    # Proportional to the square of the Fermi wavevectors
    kFp2 = (eF**2 + u**2/4) + np.sqrt(radical)
    kFm2 = (eF**2 + u**2/4) - np.sqrt(radical)

    # For kFm2, if eF > u/2, set to zero
    kFm2 = (abs(eF) <= abs(u/2)) * kFm2
    
    # Calculate the proportionality factor
    # Includes:
    #     1/(hbar vF)**2 from formula for kF
    #     1/pi from n = (kFp2 - kFm2)/pi
    #     Sets to zero if Fermi in the gap
    prop = (abs(eF)>emin(u))*np.sign(eF)*(1 / (hbar**2 * vF**2 * pi))

    return prop * (kFp2 - kFm2)




#########################
### Density of States ###
#########################
def DOS(e, u):
    '''Returns the density of states per unit area as a function of energy'''
    # DOS is symmetric for electrons and holes
    e = np.atleast_1d(abs(e))
    
    # Define the multiplicative factor out front
    # Set to 0 is energy is below the minimum
    mult = (e>emin(u)) * ( e / (pi * hbar**2 * vF**2) )
    
    # Calculate the discriminant
    # Remember, we wil need to divide by it's sqrt
    # So set values disc<=0 to 1
    # We will multiply the result by zero for these energies anyway later on.
    disc = e**2 * (g1**2 + u**2) - g1**2 * u**2 / 4
    disc = (e>emin(u))*disc + (e<=emin(u))*1
    
    # Calculate quantities proportional to derivatives of k^2
    propdkp2 = 2 + (g1**2 + u**2)/np.sqrt(disc)
    propdkm2 = 2 - (g1**2 + u**2)/np.sqrt(disc)
    
    # If energy is above sombrero region, add the positive solution
    # If within, also subtract the negative solution
    propdos = (e>emin(u))*propdkp2 - (e<=u/2)*propdkm2
    return (mult * propdos)

#####################################
### Wavefunction weights by layer ###
#####################################

def pTop_Posun(k,u):
    '''Returns the probability of finding an electron in the top layer of bilayer graphene assuming positive u and n.
    We require carrier density n because for holes, we set energy to -energy.'''
    
    # Allows us to the write rest of function assuming k is an array of values
    k = np.atleast_1d(k)
    
    # If no potential difference, return 1/2. Don't go through the song and dance.
    if u == 0:
        return 1/2 * np.ones(k.shape)
    
    # Calculate the energy
    e = Dispersion(k,u,1)
    
    # Redefine k so we don't have to write the factors all the time
    K = hbar * vF * k
    
    # Compute energy factors. pp refers to plus-plus 
    enPP = (e + u/2)**2 + K**2
    enMM = (e - u/2)**2 - K**2
    enPM = (e + u/2)**2 - K**2
    enMP = (e - u/2)**2 + K**2
    
    # Compute each layer density
    x2y2 = enPP * enMM**2 * enPM**2 / g1**2
    z2w2 = (e + u/2)**2 * enMP * enPM**2
    
    # For k=0, set x2y2 = 0 and z2w2 = 1
    x2y2 = x2y2 - (k==0)*x2y2
    z2w2 = z2w2 - (k==0)*(-1+z2w2)
    return (x2y2) / (x2y2 + z2w2)


def Pdiff(k,vminus):
    '''Returns the probability difference between finding an ELECTRON on the TOP layer minus the BOTTOM layer.'''
    
    u = -2*q*(vminus+0.0000001)
    e = Dispersion(k,u,1)
    K = hbar*vF*(k+1)
    
    numerator = (e**2 - u**2/4)**2 + 4*K**2*e**2 - K**4
    denominator = (e**2 - u**2/4)**2 + K**2*u**2 - K**4
    
    return - ( u / (2*e) ) * ( numerator / denominator )

def pTop(k,u,n):
    '''Returns the probability of finding a carrier in the top layer, regardless of whether u or n is positive or negative.
    Returns probability of finding electron in top for n>0.
    Returns probability of find hole in top for n<0.'''
    
    # Electrons are carriers (n>0), and u>0
    if n>=0 and u>=0:
        return pTop_Posun(k,u)

    # Electrons are carriers, and u<0
    # then role of top and bottom layers flip    
    if n>=0 and u<0:
        return 1 - pTop_Posun(k,abs(u))

    # Holes are carriers (n<0) and u>0
    # Prob of finding hole in top = prob of not finding an electron in the top
    if n<0 and u>=0:
        return 1 - pTop_Posun(k,u)
    
    # Holes are carriers (n<0) and u<0
    # Prob of finding hole in top, but the layers flip
    if n<0 and u<0:
        return pTop_Posun(k,abs(u))

def pBottom(k,u,n):
    return 1-pTop(k,u,n)


#########################
### Carrier Densities ###
#########################



def nplus(vplus,vminus, T, points = 10000):
    '''Returns the electron carrier density for various electrostatic potentials vplus, vminus.
    Convention is that electrons have positive carrier density while holes have negative.'''

    # Treat inputs as ndarrays so we can take advantage of broadcasting
    vplus = np.atleast_1d(vplus)
    vminus = np.atleast_1d(vminus)

    vplus = vplus.reshape(1,1,len(vplus))
    vminus = vminus.reshape(1,len(vminus),1)

    # Domain over first Brillouin zone
    ks = np.linspace(0,1/(np.sqrt(3)*a), num=points).reshape((points,1,1))

    # Calculate the kinetic energy
    KE = Dispersion(ks, -2*q*vminus,1)

    # Evaluate Fermi-Dirac
    FD = (FermiDirac(KE-q*vplus,T)-FermiDirac(KE+q*vplus,T))

    # Define integrand
    integrand = ( 2 / np.pi ) * ks * FD

    return np.squeeze(integrate.trapz(integrand,ks,axis=0))

def nminus(vplus,vminus, T,points=10000):
    '''Returns the electron carrier density for various electrostatic potentials vplus.
    Convention is that electrons have positive carrier density while holes have negative.'''

    # Treat inputs as ndarrays so we can take advantage of broadcasting
    vplus = np.atleast_1d(vplus)
    vminus = np.atleast_1d(vminus)

    vplus = vplus.reshape(1,1,len(vplus))
    vminus = vminus.reshape(1,len(vminus),1)

    # Domain over first Brillouin zone
    ks = np.linspace(0,1/(np.sqrt(3)*a), num=points).reshape((points,1,1))

    # Calculate the kinetic energy
    KE = Dispersion(ks, -2*q*vminus,1)

    # Evaluate Fermi-Dirac
    # Minus sign comes from...
    FD = (FermiDirac(KE-q*vplus,T)-FermiDirac(-KE-q*vplus,T))


    # Define integrand
    integrand =  ( 2 / np.pi ) * ks * Pdiff(ks,vminus) * FD

    return np.squeeze(integrate.trapz(integrand,ks,axis=0))