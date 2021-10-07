"""
Simple water flow example using ANUGA.
Water flowing down a channel with infiltration and changing boundary conditions
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import numpy as np
import time

import anuga
from anuga.operators.rate_operators import Rate_operator

from kostiakov import Kostiakov

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
t0 = time.time()
length = 10.
width = 5.
dx = dy = 1.           # Resolution: Length of subdivisions on both axes

domain = anuga.rectangular_cross_domain(int(length/dx), int(width/dy),
                                        len1=length, len2=width)
                                        

domain.set_name('channel2')                 # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/10                             # linear bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')  # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])   # Inflow
Br = anuga.Reflective_boundary(domain)       # Solid reflective wall
Bo = anuga.Dirichlet_boundary([-5, 0, 0])    # Outflow

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Set up infiltration
# Using the empirical Kostiakov equation Z = k.T^a + f0.T + C
# where a (dimensionless) and k (mm/h^a) are  empirical coefficients, 
#       f0 (mm/h) is a final infiltration rate
# and   C (mm) is an instantaneous adsorption and crack fill term
#------------------------------------------------------------------------------
test_infiltration = True
kost_lewis = {'a': 0.5, 'k': 40, 'f0': 0.05, 'C': 10}
if sum(kost_lewis.values()) < 0.000001:
    print("\nNo infiltration parameters - no infiltration.\n")
else:
    kost = Kostiakov(domain, kost_lewis, test_infiltration)

    def infilt_function(x, y):
        kost.mk_infilt(domain, test_infiltration)
        return -kost.inf_rate

    infilt = Rate_operator(domain, rate=infilt_function)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=40.0):
    domain.print_timestepping_statistics()
    if domain.get_quantity('stage').\
           get_values(interpolation_points=[[10, 2.5]]) > 0:        
        print('Stage > 0: Changing to outflow boundary')
        domain.set_boundary({'right': Bo})

if test_infiltration:
    kost.save_infiltration_test('') 

etxt = " Finished at {} ".format(time.strftime("%d-%b-%Y %H:%M"))
etxt += "after {:.2f} minutes...".format((time.time()-t0)/60)
print(etxt)
