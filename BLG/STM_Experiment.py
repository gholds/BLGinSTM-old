from BandStructure import *

def y(a,b,x0):
    '''Performs the integral y as defined above.
    Requires |a|>=1 and sign(a)=sign(b) with |b|>|a|.'''
    integrand = lambda t: np.cosh(t)*np.exp(np.sign(t)*np.cosh(t)/x0)
    
    limita = np.sign(a)*np.arccosh(abs(a))
    limitb = np.sign(b)*np.arccosh(abs(b))
    
    # quad returns a tuple (value,error). We only want the value.
    I = abs(integrate.quad(integrand,limita,limitb)[0])

    return I

def Y(ea,eb,u,VT):
    '''Performs the integral Y as defined above.
    Requires |ea|>emin with sign(ea)=sign(eb) and |eb|>|ea|.'''
    prefactor = emin(u)*np.sqrt(g1**2+u**2)
    
    phibar = Wtip - (q/2)*VT
    kappa0 = np.sqrt(2*m*phibar)/hbar
    x0 = 2*phibar / ( kappa0 * d1 * emin(u) )
    
    a = ea/emin(u)
    b = eb/emin(u)
    
    return prefactor * y(a,b,x0)

def Integral(ea,eb,u,VT,d):
    '''Returns the form of the integral in the mexican hat domain (d=1) or the outer domain (d=2).
    Does not include factor of 1/pi*(hbar*vF)^2.'''
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

def tunnelcurrent(eF,u,VT,s=False,WeightEnergies=False):
    '''Returns tunnel current.'''
    
    # Estimate prefactor C
    C = (4*pi*q / hbar) * 1 * Ac *1
    
    if WeightEnergies == False:
        return -C*(carrierDensity(eF-q*VT,u)-carrierDensity(eF,u))
    
    if WeightEnergies == True:
        # Calculate the parameters we need
        phibar = Wtip - (q/2)*VT
        kappa0 = np.sqrt(2*m*phibar)/hbar
        
        # Extra prefactor that came from the variable change
        C = C* np.exp(kappa0*d1*(-eF+q*VT)/(2*phibar))
        
        # Calculate the domain of integration
        # Integrate in a positive direction, then change the sign later if needed
        ea = min(eF,eF-q*VT)
        eb = max(eF,eF-q*VT)
        integrals = [0,0,0,0]

        # If ea exists in lowest domain, integrate to eb or edge of domain
        if ea < -abs(u)/2:
            integrals[0] = Integral(ea,min(eb,-abs(u)/2),u,VT,2)
            ea = -abs(u)/2
        if ea<eb and ea<-emin(u):
            integrals[1] = Integral(ea,min(eb,-emin(u)),u,VT,1)
            ea = emin(u)
        if ea<emin(u): ea=emin(u)
        if ea < eb and ea<abs(u)/2:
            integrals[2] = Integral(ea,min(eb,abs(u)/2),u,VT,1)
            ea = abs(u)/2
        if ea < eb:
            integrals[3] = Integral(ea,eb,u,VT,2)

        
        return np.sign(VT) * np.sum(integrals) / (pi*hbar**2*vF**2)

def Imeas(VT,VG,s=False,WeightEnergies=False):
    '''Returns the measured tunneling current (up to a multiplicative factor) as a function of VG and VT.
    Assumes no screening by default.'''
    
    # Get the Fermi level and interlayer potential energies
    eF = eFermi_eq(VT,VG,s=s)
    u  = u_eq(VT,VG,s=s)
    
    return tunnelcurrent(eF,u,VT,s=s,WeightEnergies=WeightEnergies)