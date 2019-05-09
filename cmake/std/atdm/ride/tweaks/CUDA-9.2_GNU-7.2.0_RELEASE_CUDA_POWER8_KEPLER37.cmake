# Disable test that times out at 10 minutes (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_main_driver_energy-ss-loca-eigenvalue_DISABLE ON)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedCurlLaplacianExample-ConvTest-Tri-Order-1_DISABLE ON)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/CUDA_COMMON_TWEAKS.cmake")

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
