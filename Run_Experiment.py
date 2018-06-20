import Data_Manipulation as dm
import numpy as np
from BLG.Universal_Constants import *
from BLG.BLG_Constants import *
from BLG.BandStructure import *
from scipy import integrate
import os


CBLG = e0 / a0

def setup_experiment(d1_set,d2_set, e1_set, e2_set, T_set, Wtip_set):
    '''Sets the variables of tunneling experiment:
    d1 : Tip-Sample distance in nm
    d2 : Backgate-Sample distance in nm
    d1 : Relative permittivity between tip and sample (vac is 1)
    d2 : Relative permittivity between backgate and sample (e.g. 3.9)
    T  : Temperature in K
    Wtip: Work function of tip in eV'''

    global d1, d2, e1, e2, T, C1, C2, Wtip

    d1 = d1_set*10**-9 # Tip-sample distance
    d2 = d2_set * 10**-9 # Backgate-Sample distance


    #e0 = 8.85 * 10**-12 # permittivity of free space
    e1 = e1_set*e0 # Vacuum between tip and sample
    e2 = e2_set*e0 # SiO2/hBn

    # Capacitances
    C1 = e1 / d1
    C2 = e2 / d2


    T=T_set

    Wtip = Wtip_set*eVtoJ



def planeplus(vplus, vminus,VGplus,VGminus):
    return (C1 + C2)*(vplus - VGplus) + (C1 - C2)*(vminus - VGminus)

def planeminus(vplus,vminus,VGplus,VGminus):
    return (C1+C2)*( vminus - VGminus ) + (C1 - C2)*( vplus - VGplus) - 4*CBLG*vminus


def generate_vplus_vminus(VTrange,num_vts_100,VBrange,num_vbs_100):
    '''Finds equilibrium values of V+ and V- for a grid of
    tip voltages VT and gate voltages VG.
    Broadcasts in batches of 10.'''

    # Load values of vplus and vminus to search over

    os.system('cls')
    print('Loading Carrier Densities')
    vplus = dm.get_vplus(T)
    vplus = vplus.reshape(1,len(vplus),1,1).astype('float32')
    vminus = dm.get_vminus(T)
    vminus = vminus.reshape(len(vminus),1,1,1).astype('float32')

    # b refers to 'batch size'
    # We loop over the range of tip and gate voltages
    # and broadcast in batches to avoid memory errors
    b = 5

    # load the carrier densities from their files
    nplus_array = dm.get_nplus(T).astype('float32')
    nminus_array = dm.get_nminus(T).astype('float32')

    nplus_array = nplus_array[:,:,np.newaxis,np.newaxis]
    nminus_array = nminus_array[:,:,np.newaxis,np.newaxis]

    # Choose tip and backgate voltages
    num_vts = int(100 * num_vts_100) # number of points for VT
    num_vbs = int(100 * num_vbs_100) # number of points for VB

    # Create array of Tip and Backgate voltages
    VT = np.linspace(VTrange[0],VTrange[1],num=num_vts)
    VT = VT.reshape(1,1,num_vts,1).astype('float32')
    VB = np.linspace(VBrange[0],VBrange[1],num=num_vbs)
    VB = VB.reshape(1,1,1,num_vbs).astype('float32')

    # Arrays where we will load the solutions
    vplus0 = np.zeros((num_vts,num_vbs))
    vminus0 = np.zeros((num_vts,num_vbs))

    print('Computing equilibrium voltages')
    for i in range(int(num_vts/b)):
        os.system('cls')
        print('{:.2}\% finished'.format( 100* i / int(num_vts/b)))
        for j in range(int(num_vbs/b)):
            
            VGplus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]+VB[:,:,:,b*j:b*j+b])
            VGminus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]-VB[:,:,:,b*j:b*j+b])

            # Generate the intersecting planes for each pair of voltages
            plane_minus_array = planeminus(vplus, vminus,VGplus,VGminus)
            plane_plus_array  = planeplus(vplus,vminus,VGplus,VGminus)

            ### GOOD UNTIL HERE ###
            # Generate the electrostatic equations
            f1 = plane_plus_array - (-q)*nplus_array
            f2 = plane_minus_array - (-q)*nminus_array

            # Find the magnitudes
            F = f1**2+f2**2

            # Minimum for each pair of voltages
            Fmins = np.min(F,axis=(0,1))

            # Array where we will put the indices of the voltages
            Fmins_args= np.empty(np.shape(Fmins),dtype=(int,2))

            # Find the indices of V+ and V-
            for k in range(b):
                for l in range(b):
                    Fmins_args[k,l] = np.where(F[:,:,k,l]==Fmins[k,l])

            # Record the values of V+ and V- based on the indices
            vplus0[b*i:b*i+b,b*j:b*j+b] = vplus[:,Fmins_args[:,:,1].flatten('C')].squeeze().reshape(b,b)
            vminus0[b*i:b*i+b,b*j:b*j+b] = vminus[Fmins_args[:,:,0].flatten('C')].squeeze().reshape(b,b)

    return (vplus0, vminus0)


#########################################
### Code for Computing Tunnel Current ###
#########################################

def y(a,b,x0):
    integrand = lambda t: np.cosh(t)*np.exp(np.sign(t)*np.cosh(t)/x0)
    
    limita = np.sign(a)*np.arccosh(abs(a))
    limitb = np.sign(b)*np.arccosh(abs(b))

    # quad returns a tuple (value,error). We only want the value.
    I = abs(integrate.quad(integrand,limita,limitb)[0])

    return I

def Y(ea,eb,u,VT):
    emin0 = emin(u)
    prefactor = emin0*np.sqrt(g1**2+u**2)
    
    phibar = Wtip - (q/2)*VT
    kappa0 = np.sqrt(2*m*phibar)/hbar
    x0 = 2*phibar / ( kappa0 * d1 * emin0 )

    a = ea/emin0
    b = eb/emin0

    # Can't call y with arguments < 1
    if abs(a) < 1 and abs(b)>=1:
        a = np.sign(b) * 1
    if abs(b) < 1 and abs(a)>=1:
        b = np.sign(a) * 1
    if a == b:
        return 0

    return prefactor * y(a,b,x0)

def Integral(ea,eb,u,VT,d):
    '''Returns the form of the integral in the mexican hat domain (d=1) or the outer domain (d=2).
    Does not include factor of 1/pi*(hbar*vF)^2.'''
    if u == 0:
        return 
    if d==1:
        return 2*Y(ea,eb,u,VT)
    if d==2:
        # Calculate the parameters we need
        phibar = Wtip - (q/2)*VT
        kappa0 = np.sqrt(2*m*phibar)/hbar
        
        integrand = lambda x: 2*abs(x)*np.exp(x*kappa0*d1/(2*phibar))
        
        return integrate.quad(integrand,ea,eb)[0]+Y(ea,eb,u,VT)
    else:
        return 0

def tunnelcurrent_obsolete(vplus,vminus,VT):
    '''Returns tunnel current.'''
    if VT==0: return 0
    eF = q*vplus
    u = -2*q*vminus

    # Estimate prefactor C
    C = (4*pi*q / hbar) * 1 * Ac *1

    # Calculate the parameters we need
    phibar = Wtip - (q/2)*VT

    kappa0 = np.sqrt(2*m*phibar)/hbar
    
    # Extra prefactor that came from the variable change
    C = C* np.exp(kappa0*d1*(-eF+q*VT)/(2*phibar))
    
    # Calculate the domain of integration
    # Integrate in a positive direction, then change the sign later if needed
    ea = min(eF,eF-q*VT)
    eb = max(eF,eF-q*VT)

    # If u==0, DOS is a constant
    if u == 0:
        dos = 2 * g1 / (pi * (hbar*vF)**2 )
        integrand = lambda x: np.exp(kappa0 * x * d1 / (2*phibar))
        return dos * integrate.quad(integrand, ea, eb)[0]
    integrals = [0,0,0,0]


    # If ea exists in lowest domain, integrate to eb or edge of domain
    if ea < -abs(u)/2:
        integrals[0] = Integral(ea,min(eb,-abs(u)/2),u,VT,2)
        ea = -abs(u)/2
    if ea<eb and ea<-emin(u):
        integrals[1] = Integral(ea,min(eb,-emin(u)),u,VT,1)
        ea = emin(u)
    if ea<emin(u):
        ea=emin(u)
    if ea < eb and ea<abs(u)/2:
        integrals[2] = Integral(ea,min(eb,abs(u)/2),u,VT,1)
        ea = abs(u)/2
    if ea < eb:
        integrals[3] = Integral(ea,eb,u,VT,2)

    
    return np.sign(VT) * np.sum(integrals) / (pi*hbar**2*vF**2)


def tunnelcurrent(vplus,vminus,VT):
    '''Returns tunnel current.'''
    if VT==0: return 0

    eF = q*vplus
    u = -2*q*vminus

    # Estimate prefactor C
    C = (4*pi*q / hbar) * 1 * Ac *1

    # Calculate the parameters we need
    phibar = Wtip - (q/2)*VT

    kappa0 = np.sqrt(2*m*phibar)/hbar
    
    # Extra prefactor that came from the variable change
    C = C* np.exp(kappa0*d1*(-eF+q*VT)/(2*phibar))
    
    # Calculate the domain of integration
    # Integrate in a positive direction, then change the sign later if needed
    ea = min(eF,eF-q*VT)
    eb = max(eF,eF-q*VT)

    integrand = lambda x : DOS(x,u) * np.exp(x*kappa0*d1/(2*phibar))

    # Points which are divergences or discontinuities
    points = np.array([u/2, -u/2, emin(u), -emin(u)])

    points = points[(ea<points) & (points<eb)]
    sign = np.sign(VT)
    return np.sign(VT) * C * integrate.quad(integrand,ea,eb,points=points)[0]
    ( pi * hbar**2 * vF**2 )


def generate_tunnelcurrent(VTrange,num_vts_100,VBrange,num_vbs_100):
    '''Generates the tunnel current over range of VTrange, VBrange 
    (lists of length 2 [min,max]. Number of point * 100.'''

    # First get equilibrium values of voltages
    vplus0, vminus0 = generate_vplus_vminus(VTrange,num_vts_100,VBrange,num_vbs_100)

    # Then generate an array of the tip voltages
    VT = np.linspace(VTrange[0],VTrange[1],num=int(100*num_vts_100))
    tc = np.empty(np.shape(vplus0))

    print('Computing tunnelcurrents')
    for i in range(np.shape(tc)[0]):
        os.system('cls')
        print('{:.2}\% finished'.format( 100* i / int(np.shape(tc)[0])))
        for j in range(np.shape(tc)[1]):
            tc[i,j] = tunnelcurrent(vplus0[i,j],vminus0[i,j],VT[i])

    return tc