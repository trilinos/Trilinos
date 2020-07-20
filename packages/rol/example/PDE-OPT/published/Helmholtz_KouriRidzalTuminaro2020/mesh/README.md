Meshes with resolutions 0 and 1 are provided as txt files.
To generate additional meshes, the meshing tool Cubit (Sandia
National Labs) and the tool ncdump must be used.

Follow these steps:

1. Modify the Nrefine parameter in the provided helmholtz.jou file,
   with 0, 1, 2, 3, 4, etc.  The provided txt mesh files correspond
   to Nrefine = 0 and Nrefine = 1.

2. Run helmholtz.jou in Cubit, and save the resulting mesh as the
   helmholtz-mesh-N.e Exodus file, where N = Nrefine, in the
   mesh subdirectory.

3. Run
   >> ncdump helmholtz-mesh-N.e > helmholtz-mesh-N.txt
