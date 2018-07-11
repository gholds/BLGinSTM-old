"""
A simple example of an animated plot
source: https://matplotlib.org/examples/animation/simple_anim.html
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from TunnelingExperiment import TunnelingExperiments
from Materials import Graphene
from UniversalConstants import q

# Set up the experiment
d1, d2 = 1, 305 # nanometers
e1, e2 = 1, 3.9 # relative permittivity
T = 0			# K
Wtip = 5 		# eV

BLG = Graphene.Bilayer()
STM = TunnelingExperiments.BLGinSTM(d1,d2,e1,e2,T,Wtip)

# Choose voltages to view
VT = np.linspace(-0.2,0.2,num=250)
VB = 30

#STM.plot_v_eq(VT,VB)


fig, (ax1,ax2,ax3,ax4) = plt.subplots( ncols=4,nrows=1,
                                       figsize=(9,7),
                                       sharey='all',
                                       gridspec_kw={'width_ratios':[1,2,2,1]})

for axis in fig.axes:
    axis.axis('off')

### Generate static dIdV at the chosen VB ###

I = np.empty_like(VT)
for i, vt in enumerate(VT):
    vp, vm = STM.v_eq(vt,VB)
    I[i] = STM.tunnelcurrent(vp,vm,vt,T)
dIdV = np.gradient(I)
np.save('dIdV_debug.npy',dIdV)


# Set the static dIdV plot
dIdV = np.load('dIdV_debug.npy')
dIdV_ax, = ax1.plot(-dIdV, VT)
ax1.set_ylim((VT[0],VT[-1]))

# Choose domain in k-space
kMax = BLG.kFermi(1*10**17, -2*q*0.1,1)
k = np.linspace(-kMax,kMax)

######
vplus, vminus = STM.v_eq(VT[0],VB)
eF, u = q*vplus, -2*q*vminus

blg_disp = BLG.Dispersion(k,u,1)/q
con = blg_disp-vplus
val = -blg_disp-vplus

dos = BLG.DOS(q*VT + eF,u)
fd  = STM.fermidirac(q*VT,VT[0])

VT_ax, = ax2.plot(k, VT[0]*np.ones_like(k))

# Plot dispersion relation
con_ax, = ax3.plot(k, con, color='k') # conduction band
val_ax, = ax3.plot(k, val, color='k') # valence band
fermi,  = ax3.plot(k, np.zeros_like(k),color='k')

dos_ax, = ax4.plot(dos/np.mean(dos), VT)
fd_ax, = ax4.plot(fd,VT)



def animate(vt):
	# Get equilibrium voltages
    vplus, vminus = STM.v_eq(vt,VB)
    eF, u = q*vplus, -2*q*vminus

    blg_disp = BLG.Dispersion(k,u,1)/q
    con = blg_disp-vplus
    val = -blg_disp-vplus

    dos = BLG.DOS(q*VT + eF,u)
    fd = STM.fermidirac(q*VT,vt)

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