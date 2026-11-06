# AMD ROCm/HIP acceleration settings for Trilinos
# Used by config/rocm/trilinos-build.cmake

# Kokkos HIP backend: enables HIP/ROCm execution space
SET_BOOL_CACHE_VAR(Kokkos_ENABLE_HIP ON)

# GPU Architecture: AMD GFX942 (MI300A APU)
# Other valid options: GFX908 (MI100), GFX90A (MI210+), GFX942 (MI300)
SET_BOOL_CACHE_VAR(Kokkos_ARCH_AMD_GFX942_APU ON)

# ROCm TPLs: AMD ROCm Libraries
SET_BOOL_CACHE_VAR(TPL_ENABLE_ROCBLAS ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_ROCSOLVER ON)
SET_BOOL_CACHE_VAR(TPL_ENABLE_ROCSPARSE ON)

SET_BOOL_CACHE_VAR(Tpetra_INST_HIP ON)
