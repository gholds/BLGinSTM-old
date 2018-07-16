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

    screening:  Boolean; "False" to neglect screening between BLG layers and "True"
                to include the effects. Calculations take much longer if "True".

    """

    def __init__(self, d1, d2, e1, e2, T, Wtip, screening=False):
        self.d1 = d1 * 10**-9
        self.d2 = d2 * 10**-9
        self.e1 = e1 * e0
        self.e2 = e2 * e0
        self.T  = T
        self.Wtip = Wtip * eVtoJ

        self.C1 = self.e1 / self.d1
        self.C2 = self.e2 / self.d2

        self.BLG = Graphene.Bilayer()

        if screening == False:
            self.screen_func = lambda x, y: y

        if screening == True:
            self.screen_func = self.BLG.screened_vminus

    def vplus_n0(self,VT,VB):
        """
        Potential of the BLG layer when no charge has accumulated

        Parameters
        ----------

        VT:     scalar or array-like, tip voltage

        VB:     scalar or array-like, gate voltage

        """
        num = self.e1*self.d2*VT + self.e2*self.d1*VB
        den = self.d1*self.e2 + self.d2*self.e1
        return num / den

    def vminus_n0(self,VT,VB):
        """
        Potential difference between the layers when no charge has accumulated.

        Parameters
        ----------

        VT:     scalar or array-like; tip voltage

        VB:     scalar or array-like; backgate voltage
        """
        num = self.BLG.d * (self.e1 + self.e2)
        den = 4*(self.e2*self.d1 + self.e1*self.d2)
        return (num/den)*(VT-VB)

    def n_exists(self,VT,VB):
        """
        Boolean function that returns whether or not charge has accumulated

        Parameters
        ---------

        VT:     scalar or array-like; tip voltage

        VB:     scalar or array-like; backgate voltage
        """

        vplus = self.vplus_n0(VT,VB)

        u = -2*q*self.vminus_n0(VT,VB)
        minimum = self.BLG.emin(u) / (q)

        return abs(vplus) >= abs(minimum)

    def nElectron(self,vplus,VT,VB):
        '''
        Returns electron density (m^-2) as a function of
        electrode-sample potential differences.

        Parameters
        ----------

        vplus:      scalar or array-like; potential of the BLG

        VT:         scalar or array-like; tip voltage

        VB:         scalar or array-like; backgate voltage

        '''
        s1 = (VT-vplus) * self.C1
        s2 = (VB-vplus) * self.C2
        
        return (s1 + s2) / q

    def vminus_n1(self,vplus,VT,VB):
        """
        The voltage difference when charge has accumulated
        """
        return (self.BLG.d / 4) * ( (VT-vplus)/self.d1 - (VB-vplus)/self.d2 )

    def vplus_root(self,vplus,VT,VB):
        """
        Self-consistency equation for the fermi level
        """
        u = -2*q*self.vminus_n1(vplus,VT,VB)
        n = self.nElectron(vplus,VT,VB)

        term1 = 4 * (q*vplus *self.BLG.g1)**2
        term2 = ( 4*(q*vplus)**2 - self.BLG.g1**2 )* u**2
        term3 = - (hbar**2 * self.BLG.vF**2 * pi)**2 * n**2

        return term1 + term2 + term3

    def vplus_n1(self,VT,VB):
        """
        Returns the Fermi level when charge has accumulated.
        Does so by finding a root to a self-consistent equation.
        """

        # otherwise we need to find the self-consistent fermi level
        f = lambda x: self.vplus_root(x,VT,VB)

        vplus = self.vplus_n0(VT,VB)

        a = min(0,vplus)
        b = max(0,vplus)

        return optimize.brentq(f,a,b)

    def v_eq(self,VT,VB):
        """
        Returns the Fermi Level first checking for charge. Uses the appropriate vplus formula.
        """

        if not self.n_exists(VT,VB):
            vp = self.vplus_n0(VT,VB)
            return (vp, self.vminus_n0(VT,VB))

        else:
            vp = self.vplus_n1(VT,VB)
            return (vp,self.vminus_n1(vp,VT,VB))

    def plot_v_eq(self,VT,VB):
        """
        Plots the equilibrium voltages of the graphene wrt the gate voltages.
        Either VT is scalar with VB array-like, or VB is scalar with VT array-like.
        """

        vp = np.empty_like(VT)
        vm = np.empty_like(VT)
        em = np.empty_like(VT)

        for i in range(len(VT)):
            vp[i], vm[i] = self.v_eq(VT[i],VB)
            em[i] = self.BLG.emin(-2*q*vm[i])/q

        plt.plot(VT,vp,'r-',label='$V_+$')
        plt.plot(VT,em,'g-',label='Emin')
        plt.plot(VT,-em,'g-')

        plt.xlabel('Tip Voltage')
        plt.ylabel('Voltages')
        plt.legend()
        plt.show()

    def generate_vplus_vminus(self,VT,VB,method):
        '''Finds equilibrium values of V+ and V- for a grid of
        tip voltages VT and gate voltages VG.

        Parameters
        ----------
        VT:     array-like, voltages over which to sweep the tip.

        VB:     array-like, voltages over which to sweep the backgate.
        
        method: Method used to compute equilibrium voltages. Only 'DasSarma' works.
        '''

        if method == 'DasSarma':

            num_vts = np.shape(VT)[0]
            num_vbs = np.shape(VB)[0]

            VT = VT.reshape((num_vts,1))
            VB = VB.reshape((1,num_vbs))

            vp = np.empty(np.shape(VT*VB))
            vm = np.empty_like(vp)

            for i in range(num_vts):
                for j in range(num_vbs):
                    if not self.n_exists(VT[i,0],VB[0,j]):
                        vp[i,j] = self.vplus_n0(VT[i,0],VB[0,j])
                        vm[i,j] = self.vminus_n0(VT[i,0],VB[0,j])
                    else:
                        vp[i,j] = self.vplus_n1(VT[i,0],VB[0,j])
                        vm[i,j] = self.vminus_n1(vp[i,j],VT[i,0],VB[0,j])

            return (vp, vm)

        if method == 'YoungLevitov':
            # Load values of vplus and vminus to search over

            print('YoungLevitov method is not supported')
            return

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

    def fermidirac(self, x, vt):
        return Temperature.FermiDirac(x-q*vt,0) - Temperature.FermiDirac(x,0)

    def tunnelcurrent(self,vplus,vminus,VT,T):
        '''Returns tunnel current.'''

        eF = -q*vplus
        u = -2*q*vminus

        # Estimate prefactor C0
        C0 = (4*pi*q / hbar) * 1 * self.BLG.Ac *1

        # Calculate the parameters we need
        phibar = self.Wtip - (q/2)*VT

        kappa0 = np.sqrt(2*m*phibar)/hbar

        fermidirac = lambda x : Temperature.FermiDirac(x-q*VT,T) - Temperature.FermiDirac(x,T)
        integrand = lambda x : fermidirac(x) * self.BLG.DOS(eF+x,u) * np.exp((x)*kappa0*self.d1/(2*phibar))

        # Points which are divergences or discontinuities or the bounds
        # At one point the end points were changing the result, but no longer
        bounds = np.sort( np.array([u/2, -u/2, self.BLG.emin(u), -self.BLG.emin(u),eF,eF+q*VT]) ) 
        bounds = bounds - eF

        tc = np.empty(len(bounds)-1)

        for i in range(len(tc)):
            tc[i] = integrate.quad(integrand,bounds[i],bounds[i+1])[0]

        return tc.sum()

    def generate_tunnelcurrent(self,VT,VB,method='DasSarma'):
        '''
        Computed tunnel current over range VT by VB.

        Parameters
        ----------
        VT:     array-like, tip voltages

        VB:     array-like, backgate voltages

        method: Method used to compute the equilibrium voltages.
                Only 'DasSarma' is valid

        
        '''

        # Get equilibrium values of voltages
        vplus0, vminus0 = self.generate_vplus_vminus(VT,VB,method)


        tc = np.empty(np.shape(vplus0))

        print('Computing tunnel currents')
        for i in range(np.shape(tc)[0]):
            print('{} % finished'.format( 100* i / int(np.shape(tc)[0])),end='\r')
            for j in range(np.shape(tc)[1]):
                tc[i,j] = self.tunnelcurrent(vplus0[i,j],vminus0[i,j],VT[i],0)

        self.VT = VT
        self.VB = VB
        self.I = tc

    def plot_dIdV(self,show=True,save=False):
        dIdV = np.gradient(self.I,axis=0) # dI/dV
        #IV = self.I / self.VT[:,np.newaxis] # I/V

        fig, ax = plt.subplots(figsize=(7,6))

        dIdV_plot = plt.imshow(dIdV,cmap=cm.RdYlGn,origin='lower',
                                aspect='auto',extent=(self.VB[0],self.VB[-1],1e3*self.VT[0],1e3*self.VT[-1]))
        fig.suptitle('$dI/dV$, Tip Height ={} nm'.format(self.d1*10**9))
        cbar = fig.colorbar(dIdV_plot,label='$dI/dV$ (S)')
        ax.set_xlabel('Gate Voltage (V)')
        ax.set_ylabel('Tip Voltage (mV)')
        if show == True:
            plt.show()

        if save == True:
            import os
            save_dir = os.path.join( os.path.dirname(__file__),
                                    'dIdV_Plots')
            fig.savefig(os.path.join(save_dir,'tip_height_{}ang.png'.format(round(self.d1*10**10))))

    def plot_schematic(self,VT, VB):
        '''
        Displays a visualization of the experiment

        Arguments
        ----------
        VT:     array-like, Range of tip voltages to sweep

        VB:     scalar, Backgate voltage at which to sweep tip voltage
        '''
        fig, (ax1,ax2,ax3,ax4) = plt.subplots( ncols=4,nrows=1,
                                               figsize=(9,7),
                                               sharey='all',
                                               gridspec_kw={'width_ratios':[1,2,2,1]})

        for axis in fig.axes:
            axis.axis('off')

        ### Generate static dIdV at the chosen VB ###

        I = np.empty_like(VT)
        for i, vt in enumerate(VT):
            vp, vm = self.v_eq(vt,VB)
            I[i] = self.tunnelcurrent(vp,vm,vt,self.T)
        dIdV = np.gradient(I)
        np.save('dIdV_debug.npy',dIdV)


        # Set the static dIdV plot
        dIdV = np.load('dIdV_debug.npy')
        dIdV_ax, = ax1.plot(-dIdV, VT)
        ax1.set_ylim((VT[0],VT[-1]))

        # Choose domain in k-space
        kMax = self.BLG.kFermi(1*10**17, -2*q*0.1,1)
        k = np.linspace(-kMax,kMax)

        vplus, vminus = self.v_eq(VT[0],VB)
        eF, u = -q*vplus, -2*q*vminus

        blg_disp = self.BLG.Dispersion(k,u,1)/q
        con = blg_disp+vplus
        val = -blg_disp+vplus

        dos = self.BLG.DOS(q*VT + eF,u)
        fd  = self.fermidirac(q*VT,VT[0])

        VT_ax, = ax2.plot(k, VT[0]*np.ones_like(k))

        # Plot dispersion relation
        con_ax, = ax3.plot(k, con, color='k') # conduction band
        val_ax, = ax3.plot(k, val, color='k') # valence band
        fermi,  = ax3.plot(k, np.zeros_like(k),color='k')

        dos_ax, = ax4.plot(dos/np.mean(dos), VT)
        fd_ax, = ax4.plot(fd,VT)

        def animate(vt):
            # Get equilibrium voltages
            vplus, vminus = self.v_eq(vt,VB)
            eF, u = -q*vplus, -2*q*vminus

            blg_disp = self.BLG.Dispersion(k,u,1)/q
            con = blg_disp+vplus
            val = -blg_disp+vplus

            dos = self.BLG.DOS(q*VT + eF,u)
            fd = self.fermidirac(q*VT,vt)

            VT_ax.set_ydata(vt*np.ones_like(k))
            ax2.fill_between(k,VT[0],vt,color='skyblue')
            con_ax.set_ydata(con)
            val_ax.set_ydata(val)

            dos_ax.set_xdata(dos/np.mean(dos))
            fd_ax.set_xdata(fd)

            return VT_ax, con_ax, val_ax, dos_ax, fd_ax,


        # Init only required for blitting to give a clean slate.
        def init():
            #dIdV.set_xdata(np.ma.array(np.sin(k), mask=True))
            VT_ax.set_ydata(np.ma.array(-VT[0]*np.ones_like(k), mask=True))
            con_ax.set_ydata(np.ma.array(con, mask=True))
            val_ax.set_ydata(np.ma.array(val, mask=True))
            dos_ax.set_xdata(np.ma.array(dos/np.mean(dos), mask=True))
            fd_ax.set_xdata(np.ma.array(fd,mask=True))
            return VT_ax, con_ax, val_ax, dos_ax, fd_ax,

        ani = animation.FuncAnimation(fig, animate, VT, init_func=init,
                                      interval=200, blit=True)
        plt.show()
        
    def planeplus(self,vplus, vminus,VGplus,VGminus):
        return (self.C1+self.C2)*(vplus - VGplus) + (self.C1-self.C2)*(vminus - VGminus)

    def planeminus(self,vplus,vminus,VGplus,VGminus):
        return (self.C1+self.C2)*(vminus - VGminus) + (self.C1 - self.C2)*(vplus - VGplus) - 4*self.BLG.C*vminus
