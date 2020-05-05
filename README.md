# Mathieu function collocation for poro-elastic plates
Code for (2D) Helmholtz scattering off multiple variable poro-elastic plates.

This folder contains code written by Matthew Colbrook for the papers:

"Fast and spectrally accurate numerical methods for perforated screens (with applications to Robin boundary conditions)"<br/>
http://www.damtp.cam.ac.uk/user/mjc249/pdfs/scattering_porous_mcolbrook.pdf<br/>
(authored with Matthew Priddin, this treated variable porosity)

and

"A Mathieu function boundary spectral method for diffraction by multiple variable poro-elastic plates, with applications to metamaterials and acoustics"<br/>
http://www.damtp.cam.ac.uk/user/mjc249/Publications.html<br/>
(authored with Anastasia Kisil, this extended to variable elastic models)

Please address questions and suggestions for additions to m.colbrook@damtp.cam.ac.uk

See license.txt.

Example.m provides an example for a single plate. Multiple plates are handled by adding to the P cell array, and Example2.m provides an example. See paper for definition of physical parameters and corresponding boundary conditions.
