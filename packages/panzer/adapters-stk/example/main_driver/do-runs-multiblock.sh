export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rcb"
#export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rib"
echo "*************************"
echo "RUNNING with 1 MPI Rank!"
echo "*************************"
mpirun -np 1 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml --kokkos-disable-warnings
echo "*************************"
echo "RUNNING with 2 MPI Ranks!"
echo "*************************"
mpirun -np 2 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml --kokkos-disable-warnings
echo "*************************"
echo "RUNNING with 3 MPI Ranks!"
echo "*************************"
mpirun -np 3 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml --kokkos-disable-warnings
echo "*************************"
echo "*************************"
