"""
Simple water flow example using ANUGA.
Water flowing down a channel with infiltration and more complex topography
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
# Setup some initial info
#------------------------------------------------------------------------------
t0 = time.time()
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/10

    N = len(x)
    for i in range(N):
        # Step
        if 10 < x[i] < 12:
            z[i] += 0.4 - 0.05*y[i]

        # Constriction
        if 27 < x[i] < 29 and y[i] > 3:
            z[i] += 2

        # Pole
        if (x[i] - 34)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 2

    return z



#------------------------------------------------------------------------------
# Setup computational domain on one processor
#------------------------------------------------------------------------------
length = 40.
width = 5.
dx = dy = .1           # Resolution: Length of subdivisions on both axes


if anuga.myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx),
                                         int(width/dy), len1=length, len2=width)
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_name('channel3')                  # Output name
    domain.set_flow_algorithm('DE0')
    domain.print_statistics()


    
    domain.set_quantity('elevation', topography)           # elevation is a function
    domain.set_quantity('friction', 0.01)                  # Constant friction
    domain.set_quantity('stage', expression='elevation')   # Dry initial condition
else:
    domain = None

#------------------------------------------------------------------------------
# Distribute domain on processor 0 to to other processors
#------------------------------------------------------------------------------
#parameters = dict(ghost_layer_width=3)
domain = anuga.distribute(domain, verbose= True)


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])          # Inflow
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Set up infiltration
# Using the empirical Kostiakov equation Z = k.T^a + f0.T + C
# where a (dimensionless) and k (mm/h^a) are  empirical coefficients, 
#       f0 (mm/h) is a final infiltration rate
# and   C (mm) is an instantaneous adsorption and crack fill term
#------------------------------------------------------------------------------
test_infiltration = True
kost_lewis = {'a': 0.03, 'k': 1.23, 'f0': 25, 'C': 20}
if sum(kost_lewis.values()) < 0.000001:
    print("\nNo infiltration parameters - no infiltration.\n")
else:
    kost = Kostiakov(domain, kost_lewis, test_infiltration)
    # if kost_lewis['k'] == 0:
    #     def infilt_function(x, y):
    #         kost.linear_infilt(domain, test_infiltration)
    #         return -kost.inf_rate
    # else:
    def infilt_function(x, y):
        kost.mk_infilt(domain, test_infiltration)
        return -kost.inf_rate

    infilt = Rate_operator(domain, rate=infilt_function)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.1, finaltime=16.0):
    if anuga.myid == 0:
        domain.print_timestepping_statistics()


    ## if domain.get_quantity('stage').\
    ##        get_values(interpolation_points=[[10, 2.5]]) > 0:
    ##     print 'Stage > 0: Changing to outflow boundary'
    ##     domain.set_boundary({'right': Bo})


domain.sww_merge(verbose=True)

anuga.finalize()


etxt = " Finished at {} ".format(time.strftime("%d-%b-%Y %H:%M"))
etxt += "after {:.2f} minutes...".format((time.time()-t0)/60)
print(etxt)
