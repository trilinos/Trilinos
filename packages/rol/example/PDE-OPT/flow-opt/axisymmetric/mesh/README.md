**Geometry and meshing**
The two problem geometries and the corresponding meshes from the paper

H. Antil, D. P. Kouri, D. Ridzal, D. B. Robinson, M. Salloum (2022),
_Uniform flow in axisymmetric devices through permeability optimization_.

are given in the Cubit journal files

  1) axi-spherical.jou

and

  2) axi-nonconvex.jou


**How to get from a Cubit journal file to a mesh text file for PDE-OPT**

1. To generate the Exodus (.exo or .e) file, on a Linux system run

     cubit -batch -nojournal -nographics -information on -warning on -input filename.jou

   (requires Cubit)

2. Then run

     ncdump filename.exo > filename.txt

   to make it into a text file readable by PDE-OPT.
   (requires ncdump)
