# -*- coding: utf-8 -*-
"""
Class definition for infitration in the Anuga 2D surface water flow model
"""
import numpy as np
import math
from os import path
#------------------------------------------------------------------------------
class Kostiakov:
    """
    a, k, f0    -- Modified Kostiakov (MK) parameters
    C           -- Initial soil infiltration for cracking soils  
    tstep       -- model timestep
    wdepth      -- surface water depth 
    mwd         -- miminum allowed water depth
    inf_t       -- infiltration opportunity (time since inundation began)
    inf_pot     -- potential infiltration depth at the start of each timestep
    inf_cum     -- cumulative infiltration
    inf_rate    -- infiltration rate (passed to the anuga rate operator) 
    """
    def __init__(self, domain, mk, test_inf=False):
        self.test = test_inf
        self.a = mk['a']
        self.k = mk['k'] * 0.001 / (3600.**self.a)    # convert mm/hr^a to m/s^a
        self.f0 = mk['f0'] * 0.001 / 3600.            # convert mm/hr to m/s
        self.C = mk['C'] * 0.001                      # convert mm to m

        self.elev = domain.get_quantity('elevation').get_values(location='centroids')
        self.inf_t = np.zeros(len(domain))
        self.inf_pot = np.zeros(len(domain))
        self.inf_cum = np.zeros(len(domain))
        self.inf_rate = np.zeros(len(domain))
        # avoids a second power calculation each timestep
        self.inft_pow = np.zeros(len(domain))

        self.mwd = domain.get_minimum_allowed_height()
        if test_inf:
            c_cords = domain.get_centroid_coordinates(absolute=True)
            xcoord = c_cords[:,0]
            ycoord = c_cords[:,1]
            poly_extent = np.array(domain.get_extent())
            x1= (poly_extent[1] - poly_extent[0])/2
            y1= (poly_extent[3] - poly_extent[2])/4
            short_d = 1.0e+100
            for k in range(len(domain)):
                if math.hypot(xcoord[k] - x1, ycoord[k] - y1) < short_d:
                    self.index_k = k
                    short_d = math.hypot(xcoord[k] - x1, ycoord[k] - y1)

            self.x = xcoord[self.index_k]
            self.y = ycoord[self.index_k]
            self.cell_oppt = []
            self.cell_inf = []
            self.model_time = []
            self.cell_unsat = []
            self.cell_wdepth = []

    def mk_infilt(self, domain, test=False): 
        """
        Simulates the removal of surface water from the model domain at a rate
        calculated with the modified Kostiakov-Lewis (MK) infiltration equation
        (Hartley, 1992) with an additional term C that accounts for
        "instantaneous" crack fill (Austin & Prendergast, 1997).

        The MK equation describes infiltration under the condition of abundant
        surface water, which is not the case at the advancing front of overland
        flow where infiltration rate can be constrained to the rate at which
        surface water arrives.
        
        At each model timestep and for each each cell in the domain the algorithm
        calculates the potential infiltration on "wet" cells. Each cell must be
        individually tracked because infiltration will start and end at different
        times on each cell.
        
        On cells where inf_pot is less than water depth (wdepth), the infiltration
        rate (inf_rate) is inf_pot/tstep and because inf_pot is satisfied, inf_t
        is incremented. Conversely, on cells where water depth is less than
        inf_pot, inf_rate is wdepth/tstep, inf_pot is reduced by wdepth and inf_t
        is not incremented.
        """
        tstep = domain.get_timestep()
        wdepth = domain.get_quantity('stage').get_values(location='centroids') - self.elev
        
        self.inf_rate = np.zeros(len(wdepth))

        # Initiate infiltration for cells where surface water exists and
        # infiltration opportunity time = 0,
        idx = np.where(np.logical_and(wdepth >= self.mwd, self.inf_t == 0))[0] 
        if len(idx) > 0:
            self.inf_t[idx] = tstep
            self.inft_pow[idx] = self.inf_t[idx] ** self.a
            self.inf_pot[idx] = self.k * self.inft_pow[idx] \
                    + self.f0 * self.inf_t[idx] + self.C

        # On cells where where surface water exists and the infiltration potential
        # has been met, increment the infiltration time step and calculate the new
        # infiltration potential. 
        idx1 = np.where(np.logical_and(wdepth >=self.mwd, self.inf_pot == 0))[0] 
        if len(idx1) > 0:
            inft = self.inft_pow[idx1]
            self.inf_t[idx1] += tstep
            self.inft_pow[idx1] = self.inf_t[idx1] ** self.a
            self.inf_pot[idx1] = self.k * (self.inft_pow[idx1] - inft) + self.f0 * tstep
                            
        idx2 = np.where(np.logical_and(wdepth >= self.mwd, self.inf_pot < wdepth))[0] 
        idx3 = np.where(np.logical_and(wdepth >= self.mwd, self.inf_pot >= wdepth))[0]

        # Calculate cell infiltration rate and update potential infiltration
        if len(idx2) > 0:
            self.inf_rate[idx2] = self.inf_pot[idx2] / tstep
            self.inf_pot[idx2] = 0.

        if len(idx3) > 0:
            self.inf_rate[idx3] = wdepth[idx3] / tstep
            self.inf_pot[idx3] -= wdepth[idx3]
        
        # Update cumulative infiltration
        self.inf_cum += self.inf_rate * tstep

        if test:
            self.cell_oppt.append(self.inf_t[self.index_k])
            self.cell_inf.append(self.inf_cum[self.index_k])
            self.cell_unsat.append(self.inf_pot[self.index_k])
            self.cell_wdepth.append(wdepth[self.index_k])
            self.model_time.append(domain.get_time())

    def save_infiltration_test(self, pth):
        """
        Writes infiltration output for a single mesh cell to cell_infilt.csv in the
        scenario directory.

        Parameters:
            pth: the directory path
        """
        fn = path.join(pth + 'cell_infilt' + str(self.index_k) + '.csv')
        with open(fn , 'w') as fw:
            fw.write('Model_time(s), Opp_time(s), Infiltration(m), Inf_pot(m), Wdepth(m)')
            fw.write(',,,,,, %s, %s, %s\n' % (str(self.index_k), self.x, self.y))
            for f1, f2, f3, f4, f5 in zip(self.model_time, self.cell_oppt,
                                            self.cell_inf, self.cell_unsat, self.cell_wdepth):
                fw.write('%s,%s,%s,%s,%s\n' % (f1, f2, f3, f4, f5))

#------------------------------------------------------------------------------
if __name__ == "__main__":
    pass
