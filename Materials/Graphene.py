from . import *     # Import from __init__.py
from abc import ABCMeta # For inheritance
from BLG.Universal_Constants import *

class BaseGraphene:
    """
    Base class for all types of graphene


    """

    __metaclass__ = ABCMeta

    def __init__(self):
        self.a = 1.42 * (10**(-10)) # (m), Interatom spacing
        self.Ac = 3*np.sqrt(3)*(self.a**2) / 2 # (m^2), Area of unit cell of graphene
        self.g0 = 2.8*eVtoJ # (J), Interatom hopping potential
        self.vF = 3*self.a*self.g0/(2*hbar)

    def FermiDirac(self, e, T):
        """
        Returns the Fermi-Dirac distribution.

        Parameters
        ----------
        e   :   Energy in J

        T   :   Temperature (K)
        """
        return np.exp(-np.logaddexp(e/(kB*(T+0.0001)),0))


class Bilayer(BaseGraphene):
    """
    Bilayer graphene

    """

    g1  = 0.358*eVtoJ # (J), A1-B1 hopping potential of BLG
    g3  = 0.3 * eVtoJ # (J), A1-B2 hopping potential of BLG
    d   = 3*(10**-10) # (m), interlayer spacing of BLG
    approx_choices = ['None', 'g3=0', 'LowEnergy', 'Quadratic']


    def Hamiltonian(self,k,u):
        '''
        Returns the full tight-binding Hamiltonian of BLG.
        Does not support vectorization.

        Parameters
        ----------
        k:  Wavenumber (1/m).

        u:  Interlayer potential energy difference (J).
        '''
        H = np.array([  [-u/2, hbar*self.vF*k,0,3*self.g3*self.a*k],
                [hbar*self.vF*k, -u/2, self.g1, 0   ],
                [0, self.g1, u/2, hbar*self.vF],
                [3*self.g3*self.a*k, 0, hbar*self.vF*k,0]])
        return H
    ######################
    ### Band Structure ###
    ######################

    def Dispersion(self,k,u,band,approx='g3=0'):
        '''
        Returns the energy (J) of an electron with wavevector k (rad/m)
        in first (band=1) or second (band=2) conduction band.
        Only approximation is g3=0.
        To get holes, simply result multiply by -1.
        '''
        p = hbar * self.vF * k

        if approx == 'g3=0':
            radical=(self.g1**4)/4 + (u**2 + self.g1**2)*(p**2)
            return np.sqrt( (self.g1**2)/2 + (u**2)/4 + p**2 + ((-1)**(band))*np.sqrt(radical) )

        if approx == 'LowEnergy':
            '''
            Low Energy effective. Eigenvalues of 

            H = ( ( u/2, p^2 / 2m ) , ( p^2/2m, -u/2 ) )

            '''

            meff = ( self.g1 / (2 * (self.vF)**2) )**-1
            return np.sqrt( (u/2)**2 + ( (hbar * k)**2 / 2*meff )**2 )

        if approx == 'None':
            '''
            No approximation. Compute eigenvalues of Hamiltonian
            '''

            H = np.array([  [-u/2, hbar*self.vF*k,0,3*self.g3*self.a*k],
                            [hbar*self.vF*k, -u/2, self.g1, 0   ],
                            [0, self.g1, u/2, hbar*self.vF],
                            [3*self.g3*self.a*k, 0, hbar*self.vF*k,0]])
            return np.abs(linalg.eigvalsh(H)).min()

    def kmin(u, band=1):
        '''
        Returns positive wavenumber at the minimum of the first band in 1/m.
        '''
        k2 = ( u**2 / (2*hbar*vF)**2 ) * ( (2*g1**2 + u**2) /( g1**2 + u**2 ) )
        return np.sqrt(k2)

    def emin(self,u):
        '''
        Returns minimum of the first band in Joules.
        '''
        emin2 = (u/2)**2 * (self.g1**2 / ( self.g1**2 + u**2 ) )
        return np.sqrt(emin2)

    def DOS(self, e, u):
        '''
        Returns the density of states per unit area (1/m^2) as 
        a function of energy given the gap u
        '''
        e = np.atleast_1d(abs(e))
        
        # Define the multiplicative factor out front
        # Set to 0 is energy is below the minimum
        mult = (e>self.emin(u)) * ( e / (pi * hbar**2 * self.vF**2) )
        
        # Calculate the discriminant
        # Remember, we wil need to divide by it's sqrt
        # So set values disc<=0 to 1
        # We will multiply the result by zero for these energies anyway later on.
        disc = e**2 * (self.g1**2 + u**2) - self.g1**2 * u**2 / 4
        disc = (e>self.emin(u))*disc + (e<=self.emin(u))*1
        
        # Calculate quantities proportional to derivatives of k^2
        propdkp2 = 2 + (self.g1**2 + u**2)/np.sqrt(disc)
        propdkm2 = 2 - (self.g1**2 + u**2)/np.sqrt(disc)
        
        # If energy is above sombrero region, add the positive solution
        # If within, also subtract the negative solution
        propdos = (e>self.emin(u))*propdkp2 - (e<=u/2)*propdkm2
        return (mult * propdos)

    def Pdiff(self,k,vminus,approx='g3=0'):
        '''Returns the probability difference between finding an ELECTRON on the TOP layer minus the BOTTOM layer.'''
        
        u = -2*q*(vminus+0.0000001)
        
        if approx=='g3=0':
            e = self.Dispersion(k,u,1)
            K = hbar*self.vF*(k+1)
            
            numerator = (e**2 - u**2/4)**2 + 4*K**2*e**2 - K**4
            denominator = (e**2 - u**2/4)**2 + K**2*u**2 - K**4
            
            return - ( u / (2*e) ) * ( numerator / denominator )

        if approx=='LowEnergy':
            meff = ( self.g1 / (2 * (self.vF)**2) )
            print(meff)
            denominator_squared = ( ( (hbar*k)**2/meff )**2 + u**2 )
            
            return - u / np.sqrt(denominator_squared)

    def kFermi(self,n,u,pm):
        '''
        Returns Fermi vector kF+ for pm=1 and kF- for pm=2 in units rad/m
        '''
            
        # Define the more complicated factors and terms
        numerator = (pi * hbar**2 *self.vF**2 * n)**2 + ( self.g1*u )**2
        denominator = self.g1**2 + u**2
        pmterm = 2*pi*hbar**2 * self.vF**2 * abs(n) # abs for fact that electrons and holes symmetric
        
        # Factor proportional to k**2
        propk2 = ( numerator / denominator ) + u**2 + (-1)**(pm-1) * pmterm
        
        # If the fermi level is above u/2, set kF- to zero
        # This says that the region of occupied states is now a disk
        if pm%2==0:
            propk2 = (propk2 >= 0) * propk2
            propk2 = (Dispersion(kFermi(n,u,1),u,1)<u/2) * propk2
        
        return np.sqrt( propk2 ) / (2*hbar*self.vF)

    def eFermi(self,n,u):
        '''
        Returns the Fermi level (Joules) given density n and interlayer potential energy difference u
        Positive n returns a positive Fermi level, meaning positive carrier densities are electrons by convention.
        '''
        
        numerator = (hbar**2 * self.vF**2 * n *pi)**2 + (self.g1 * u)**2
        denominator = 4 * (self.g1**2 + u**2)
        
        return np.sign(n) * np.sqrt( numerator / denominator )

    #########################
    ### Carrier Densities ###
    #########################

    def nplusT0(self,vplus,vminus,approx='g3=0'):
        """
        Analytically computes the electron density at zero temperature.
        Faster than Bilayer.nplus() since this function allows
        for vectorized operations.
        """

        # Convert voltages to energies
        eF = eVtoJ*vplus
        u  = -2*eVtoJ*vminus

        if approx == 'g3=0':
            # Calculate the radical
            radical = (self.g1**2+u**2) * eF**2 - self.g1**2 * u**2 / 4

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
            prop = (abs(eF)>self.emin(u))*np.sign(eF)*(1 / (hbar**2 * self.vF**2 * pi))

            return prop * (kFp2 - kFm2)

        if approx == 'LowEnergy':
            """
            See Young and Levitov 2011.
            """
            meff = ( self.g1 / (2 * (self.vF)**2) )

            nu0 = 2 * meff * q  / (pi * hbar**2)

            energy_diff = (np.abs(eF)>np.abs(u/2)) * (eF**2 - (u/2)**2)
            return (nu0/q) * np.sign(eF) * np.sqrt(energy_diff)

    def nplus(self,vplus,vminus, T, approx='g3=0',points = 10000):
        '''
        Returns the electron carrier density for various electrostatic potentials vplus, vminus.
        Convention is that electrons have positive carrier density while holes have negative.
        '''

        # Treat inputs as ndarrays so we can take advantage of broadcasting
        vplus = np.atleast_1d(vplus)
        vminus = np.atleast_1d(vminus)

        vplus = vplus.reshape(1,1,len(vplus))
        vminus = vminus.reshape(1,len(vminus),1)

        # Domain over first Brillouin zone
        ks = np.linspace(0,1/(np.sqrt(3)*self.a), num=points).reshape((points,1,1))

        # Calculate the kinetic energy
        KE = self.Dispersion(ks, -2*q*vminus,1,approx)

        # Evaluate Fermi-Dirac
        FD = (self.FermiDirac(KE-q*vplus,T)-self.FermiDirac(KE+q*vplus,T))

        # Define integrand
        integrand = ( 2 / np.pi ) * ks * FD

        return np.squeeze(integrate.trapz(integrand,ks,axis=0))

    def nminus(self,vplus,vminus, T, approx='g3=0', points=10000):
        '''
        Returns the electron carrier density for various electrostatic potentials vplus.
        Convention is that electrons have positive carrier density while holes have negative.
        '''
        # Treat inputs as ndarrays so we can take advantage of broadcasting
        vplus = np.atleast_1d(vplus)
        vminus = np.atleast_1d(vminus)

        vplus = vplus.reshape(1,1,len(vplus))
        vminus = vminus.reshape(1,len(vminus),1)

        # Domain over first Brillouin zone
        ks = np.linspace(0,1/(np.sqrt(3)*self.a), num=points).reshape((points,1,1))

        # Calculate the kinetic energy
        KE = self.Dispersion(ks, -2*q*vminus,1, approx)

        # Evaluate Fermi-Dirac
        # Minus sign comes from...
        FD = (self.FermiDirac(KE-q*vplus,T)-self.FermiDirac(-KE-q*vplus,T))


        # Define integrand
        integrand =  ( 2 / np.pi ) * ks * self.Pdiff(ks,vminus,approx) * FD

        return np.squeeze(integrate.trapz(integrand,ks,axis=0))
