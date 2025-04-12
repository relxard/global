# global
A finite difference solution of the  advection-diffusion equation on a global latitude-longitude grid.
The model can run on any 1 or 2.5 deg grid with a regular grid spacing. The concentration grid is 
automatically configured to match the meteorological grid.  The grid system indicies increase from
south to north and from west to east. The longitude coordinate system can go from 000 to 360 or 
from -180 to 180. The first and last longitude grid points are adjacent to each other but do not 
overlap (cyclic boundary condition). The south to north coordinate system would run from -90 
to +90 degrees latitude. The grid's vertical coordinate system is on pressure-sigma surfaces. 
Pollutants are dispersed and advected in terms of their mass mixing ratio which is converted 
to concentration at STP for output.
