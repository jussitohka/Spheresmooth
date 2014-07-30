Spheresmooth
============

Surface smoothing with sphere prior

Author Jussi Tohka, Institute of Signal Processing,
Tampere University of Technology, Finland

This open-source software implements the surface smoothing algorithm described in
J.~Tohka. Surface smoothing based on a sphere shape model.
In J.~Tanskanen, editor, Proc. Nordic Signal Processing Symposium
(NORSIG04)}, pages 17 -- 20. IEEE 2004. It is written in C.  

http://legacy.spa.aalto.fi/sig-legacy/norsig2004/publications/12_TOHKA.PDF

The
smoothing algorithm works by 
1) converting a triangle mesh into the
dual simplex mesh, 
2) smoothing the dual simplex mesh, and 
3) converting the smoothed simplex mesh back to triangular mesh. 

The triangle meshes are given in .obj format.  
See http://www.bic.mni.mcgill.ca/users/david/FAQ/polygons_format.txt
for the file format. The program accepts only ASCII .obj files and it
does not allow giving different colors to different vertices. Not much
checking is done, so crashes may follow if an incompatible file is
given as an input.

IMPORTANT: The triangle meshes MUST BE proper simplicial
complexes. Otherwise the conversion to simplex mesh won't work. This
is a feature not a bug. The
trouble is that many isosurface algorithms (e.g. Matlab's isosurface
function) have tendency to produce meshes which are not proper
simplicial complexes.  

-----------------------------------------------------------
INSTALLATION
-----------------------------------------------------------
Just unzip the zip-file into the directory you want to install and
type "make". The software has been tested with gcc (versions gcc (GCC)
3.3.1 (Mandrake Linux 9.2 3.3.1-2mdk), and gcc 4.2.1) running on
Linux.  The code should compile on other platforms and compilers
but that has not been tested.

-----------------------------------------------------------
USAGE 
-----------------------------------------------------------
 
Usage: spriorsmooth infile.obj outfile.obj [lambda] [multiplier]

The parameter lambda specifies the amount of smoothing (default 400)
The parameter multiplier specifies the desired accuracy of computation
(default 0.1), increment the value for a speed up.

-----------------------------------------------------------
NOTES
-----------------------------------------------------------
 - The convergence issue has been solved. See the file convergence.pdf
for details. The parameter "eta" is now selected optimally.
