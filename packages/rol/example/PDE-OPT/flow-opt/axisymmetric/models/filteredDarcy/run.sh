rm control.txt
rm target.txt
rm permeability*
rm velocity*
mpirun -np 20 ./ROL_example_PDE-OPT_flow-opt_axisymmetric_models_filteredDarcy_example_01.exe display
cat permeability_* > permeability.txt
cat velocity_* > velocity.txt
