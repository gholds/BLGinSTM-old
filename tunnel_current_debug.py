import numpy as np
import Run_Experiment as re
from BLG.BandStructure import emin
from BLG.BLG_Constants import eVtoJ,q,m,hbar

d1 = 1
d2 = 305
Wtip = 5

re.setup_experiment(d1, d2, 1, 3.9, 0, Wtip)

## Full Test
VTmin , VTmax = -0.2, 0.2
num_vts_100 = 0.1 # 100 points


VBmin, VBmax = -45, 45
num_vbs_100 = 0.5 # 10 points

I = re.generate_tunnelcurrent([VTmin,VTmax], num_vts_100, [VBmin,VBmax], num_vbs_100)
