# Anuga infiltration

This code simulates infiltration for the Anuga 2D hydrodynamic model.
Anuga can be downloaded from
       <https://github.com/GeoscienceAustralia/anuga_core>
Class Kostiakov simulates the removal of surface water from the Anuga
model domain at a rate calculated with the modified Kostiakov-Lewis
(MK) infiltration equation (Hartley, 1992) with an additional term C
that accounts for "instantaneous" crack fill (Austin & Prendergast, 1997).

Class Kostiakov calculates the infiltration rate for each model cell at
each model timestep.  The rate can then be passed to a standard ANUGA
rate operator.

    Infiltrated depth Z = k * T ** a + f0 * T + C

    where a and k are empirical coefficients,
          f0 approximates the final infiltration rate,
      and C is instantaneous crack fill.

Any (or all) of these parameters may be zero.
  
The MK equation describes infiltration under the condition of abundant
surface water, which is usually not the case at the advancing front of
overland flow where infiltration rate can be constrained to the rate at
which surface water arrives.

At each model timestep and for each each cell in the model domain the
algorithm calculates the potential infiltration on "wet" cells. Each
cell must be individually tracked because infiltration will start and
end at a different time on each cell.

On cells where the water depth (wdepth) is greater than the infiltration
potential (inf_pot), the infiltration rate is inf_pot/tstep and because
inf_pot is satisfied, infiltration opportunity time (inf_t) is
incremented. Conversely, on cells where water depth is less than
inf_pot, inf_rate is wdepth/tstep, inf_pot is reduced by wdepth, wdepth
becomes 0, inf_t is not incremented and because wdepth=0 surface flow
from the cell will not occur.

Class Kostiakov provides an option to test its operation by reporting
infiltration on an single cell for all model time steps.

References:

Austin, N.R., Prendergast, J.B., 1997. Use of kinematic wave theory to
 model irrigation on cracking soil. Irrigation Science 18, 1–10.

Githui, F., Hussain, A., Morris, M., 2015. Adapting ANUGA model for
 border-check irrigation simulation, in: 21st International Congress
 on Modelling and Simulation (MODSIM), 29 Nov to 4 Dec 2015, Gold
 Coast, Australia.

Githui, F. & Hussain, A. & Morris, M. 2020. Incorporating infiltration
 in the two-dimensional ANUGA model for surface irrigation simulation.
 Irrigation Science. 38.

Hartley, D.M., 1992. Interpretation of Kostiakov Infiltration Parameters
 for Borders. J. Irrig. Drain Eng. 118, 156–165.

Morris, M., Githui, F., Hussain, A., 2015. Application of Anuga as a 2D
 surface irrigation model, 21st International Congress on Modelling and
 Simulation (MODSIM), 29 Nov to 4 Dec 2015, Gold Coast, Australia.
