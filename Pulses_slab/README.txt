Description: 
The program pulses_slab calculates the temperature evolution of an
X-ray pulse train hitting a disk in normal incidence. 
The idl program pulses_slab.pro is a copy with minor changes 
of the 2007 program pulse_laue.pro. 
With small adjustments, pulses_slab_gdl runs with 
the free software gdl (and also still with idl). 
The code pulses_slab.m is a translation into Matlab and was 
developed and tested on the free software Octave. 
The speed of the idl/gdl codes can be improved, if the 
range of the spline for the temperature axis is adjusted 
to the simulation range and therefore the number of spline 
interpolation points is reduced.   
The original code pulses_slab.pro is the most accurate, because 
the heat conduction and the specific heat are calculated at edges and 
centers of the volume elements, respectively. 
In order to speed up the codes, the Matlab and GDL codes calculate 
the ratio in the middle of each cell, which leads to less precise 
conservation of the total energy. 

benchmarks
10+2 pulses, start at 100K		time	
pulses_slab.pro on idl/Palx		13 sec  (second run, ends at 173.042K)
pulses_slab.pro on idl/exflserv02	31 sec (second run, ends at 173.042 K)
pulses_slab.m on Matlab/exflserv02	87 sec (173.043K)
pulses_slab.pro on gld/Mac		900 sec (173.042K)
pulse_slab.m on Octave/Mac		900 sec (173.043K)

