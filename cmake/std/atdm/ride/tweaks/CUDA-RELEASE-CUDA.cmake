INCLUDE("${CMAKE_CURRENT_LIST_DIR}/CUDA_COMMON_TWEAKS.cmake")

# Disable test that times out at 10 minutes (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_main_driver_energy-ss-loca-eigenvalue_DISABLE ON)

# This test fails consistently with a major numerical error (#2474)
ATDM_SET_ENABLE(Piro_MatrixFreeDecorator_UnitTests_MPI_4_DISABLE ON)
