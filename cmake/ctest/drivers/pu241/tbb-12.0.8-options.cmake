
SET(INTEL_12_0_8_ROOT /opt/intel-12.0.8)
SET(TBBROOT ${INTEL_12_0_8_ROOT}/tbb)

SET(TBB_INCLUDE_DIRS  ${TBBROOT}/include                                      CACHE FILEPATH "")
SET(TBB_LIBRARY_DIRS  ${TBBROOT}/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21  CACHE FILEPATH "")
