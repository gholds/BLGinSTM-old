from . import *
from Materials import Graphene

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
                print('{} % finished'.format( 100*(i_frac+j_frac) ) )
                VGplus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]+VB[:,:,:,b*j:b*j+b])
                VGminus = (1 / 2) * (VT[:,:,b*i:b*i+b,:]-VB[:,:,:,b*j:b*j+b])

                # Generate the intersecting planes for each pair of voltages
                plane_minus_array = self.planeminus(vplus, vminus,VGplus,VGminus)
                plane_plus_array  = self.planeplus(vplus,vminus,VGplus,VGminus)

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
