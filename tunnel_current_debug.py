import numpy as np
import Run_Experiment as re

re.setup_experiment(1, 305, 1, 3.9, 0, 5)

VTmin , VTmax = -0.2, 0.2
num_vts_100 = 0.05 # 100 points


VBmin, VBmax = -45, 45
num_vbs_100 = 1 # 100 points

I = re.generate_vplus_vminus([VTmin,VTmax], num_vts_100, [VBmin,VBmax], num_vbs_100)
