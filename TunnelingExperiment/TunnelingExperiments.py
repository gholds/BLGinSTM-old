from . import * # Packages imported in __init__.py

class BLGinSTM:

    """
    A Tunneling experiment containing information to calculate tunnel
    current from BLG.

    Parameters
    ----------

    d1, d2  :   The top gate to sample and bottom gate to sample distances
                respectively in nanometers. Converted to meters for calcs.
    
    e1, e2  :   The relative permittivities in the top gate-sample and
                bottom gate-sample regions respectively. Multiplied by e0
                for calculations.

    T       :   Temperature (K)

    Wtip    :   Work function of top gate in eV. Converted to J for calculations

    material:   A material chosen from the Materials package. So far, only bilayer
                graphene is available.

    """

    def __init__(self, d1, d2, e1, e2, T, Wtip):
        self.d1 = d1 * 10**-9
        self.d2 = d2 * 10**-9
        self.e1 = e1 * e0
        self.e2 = e2 * e0
        self.T  = T
        self.Wtip = Wtip * eVtoJ

        self.C1 = self.e1 / self.d1
        self.C2 = self.e2 / self.d2

        self.BLG = Graphene.Bilayer()


    def planeplus(self,vplus, vminus,VGplus,VGminus):
        return (self.C1+self.C2)*(vplus - VGplus) + (self.C1-self.C2)*(vminus - VGminus)

    def planeminus(self,vplus,vminus,VGplus,VGminus):
        return (self.C1+self.C2)*(vminus - VGminus) + (self.C1 - self.C2)*(vplus - VGplus) - 4*self.BLG.C*vminus


    def generate_vplus_vminus(self,VTrange,num_vts_100,VBrange,num_vbs_100):
        '''Finds equilibrium values of V+ and V- for a grid of
        tip voltages VT and gate voltages VG.
        Broadcasts in batches of 10.'''

        # Load values of vplus and vminus to search over

        print('Loading Carrier Densities')
        vplus = self.BLG.get_vplus(self.T)
        vplus = vplus.reshape(1,len(vplus),1,1).astype('float32')
        vminus = self.BLG.get_vminus(self.T)
        vminus = vminus.reshape(len(vminus),1,1,1).astype('float32')

        # b refers to 'batch size'
        # We loop over the range of tip and gate voltages
        # and broadcast in batches to avoid memory errors
        b = 5

        # load the carrier densities from their files
        nplus_array = self.BLG.get_nplus(self.T).astype('float32')
        nminus_array = self.BLG.get_nminus(self.T).astype('float32')

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
            i_frac = i / int(num_vts/b)
            
            for j in range(int(num_vbs/b)):
                j_frac = ( j / int(num_vbs/b) ) * (1/int(num_vts/b))
                print('{} % finished'.format( round(100*(i_frac+j_frac)) ), end='\r' )
                VGplus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]+VB[:,:,:,b*j:b*j+b])
                VGminus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]-VB[:,:,:,b*j:b*j+b])

                # Generate the intersecting planes for each pair of voltages
                plane_minus_array = self.planeminus(vplus, vminus,VGplus,VGminus)
                plane_plus_array  = self.planeplus(vplus,vminus,VGplus,VGminus)

                # Generate the electrostatic equations
                f1 = plane_plus_array - (-q)*nplus_array
                f2 = plane_minus_array - (-q)*nminus_array

                # Find the magnitudes
                F = f1**2 + f2**2

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

        print('100 % Finished')
        return (vplus0, vminus0)


    def tunnelcurrent(self,vplus,vminus,VT):
        '''Returns tunnel current.'''
        if VT==0: return 0

        eF = q*vplus
        u = -2*q*vminus

        # Estimate prefactor C
        C = (4*pi*q / hbar) * 1 * self.BLG.Ac *1

        # Calculate the parameters we need
        phibar = self.Wtip - (q/2)*VT

        kappa0 = np.sqrt(2*m*phibar)/hbar
        
        # Extra prefactor that came from the variable change
        C = C* np.exp(kappa0*self.d1*(-eF+q*VT)/(2*phibar))
        
        # Calculate the domain of integration
        # Integrate in a positive direction, then change the sign later if needed
        ea = min(eF,eF-q*VT)
        eb = max(eF,eF-q*VT)

        integrand = lambda x : self.BLG.DOS(x,u) * np.exp(x*kappa0*self.d1/(2*phibar))

        # Points which are divergences or discontinuities
        points = np.array([u/2, -u/2, self.BLG.emin(u), -self.BLG.emin(u)])

        points = points[(ea<points) & (points<eb)]
        sign = np.sign(VT)
        return np.sign(VT) * C * integrate.quad(integrand,ea,eb,points=points)[0]
        ( pi * hbar**2 * vF**2 )

    def generate_tunnelcurrent(self,VTrange,num_vts_100,VBrange,num_vbs_100):
        '''Generates the tunnel current over range of VTrange, VBrange 
        (lists of length 2 [min,max]. Number of point * 100.'''

        # First get equilibrium values of voltages
        vplus0, vminus0 = self.generate_vplus_vminus(VTrange,num_vts_100,VBrange,num_vbs_100)

        # Then generate an array of the tip voltages
        VT = np.linspace(VTrange[0],VTrange[1],num=int(100*num_vts_100))
        tc = np.empty(np.shape(vplus0))

        print('Computing tunnelcurrents')
        for i in range(np.shape(tc)[0]):
            print('{:.2}\% finished'.format( 100* i / int(np.shape(tc)[0])))
            for j in range(np.shape(tc)[1]):
                tc[i,j] = self.tunnelcurrent(vplus0[i,j],vminus0[i,j],VT[i])

        return tc