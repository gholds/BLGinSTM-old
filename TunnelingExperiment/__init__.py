import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
#from matplotlib.patches import ConnectionPatch
from scipy import integrate, optimize
from StatisticalDistributions import Temperature
from UniversalConstants import *
from Materials import Graphene

machine_epsilon = 7./3 - 4./3 - 1