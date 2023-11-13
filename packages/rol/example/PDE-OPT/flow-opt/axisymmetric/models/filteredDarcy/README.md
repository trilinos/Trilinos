**Filtered Darcy**

This module is used to design spatially varying porosity profiles in columns
that can be additively manufactured (printed in three dimensions) and are
used for chemical reactions.  The module is accompanied by the paper

H. Antil, D. P. Kouri, D. Ridzal, D. B. Robinson, M. Salloum (2022),
_Uniform flow in axisymmetric devices through permeability optimization_.

-- To run the code, first compile it using the standard Trilinos
configure/build process.  Make sure to configure with the full Tpetra solver
stack enabled.

-- Second, relative to the directory containing this README file, go to the
directory ../../mesh and follow the instructions to generate a mesh for the
simulation and optimization.

-- Third, run

>> ./run.sh

which executes

  rm control.txt
  rm target.txt
  rm permeability*
  rm velocity*
  mpirun -np 20 ./ROL_example_PDE-OPT_flow-opt_axisymmetric_models_filteredDarcy_example_01.exe display
  cat permeability_* > permeability.txt
  cat velocity_* > velocity.txt

Change the number of processes from 20 to whatever suits your hardware.

Note that the executable reads the file input.xml, which is identical
to the file input-spherical.xml as distributed.  To run a different example,
modify input.xml.  For instance, copy input-nonconvex.xml into input.xml,
or generate an entirely different problem setup.

-- Fourth, postprocess the generated data by running the Matlab script

  plotDarcyAxi.m

which takes in a single parameter, 1 or 0, to save or not save the figures.


Note: The code is capable of continuation in the sense that repeatedly running

  mpirun -np 20 ./ROL_example_PDE-OPT_flow-opt_axisymmetric_models_filteredDarcy_example_01.exe display

will use the stored control variables, in the files control.txt and target.txt,
as the initial guesses for subsequent runs.  To start from scratch, delete the
files control.txt and target.txt.
