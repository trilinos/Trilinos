export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rcb"
#export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rib"
mpirun -np 3 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml
mpirun -np 2 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml
mpirun -np 1 ./PanzerAdaptersSTK_main_driver.exe --i=energy-ss-tp-multiblock-ic-bc-issue.xml
