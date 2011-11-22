
SET(INTEL_12_0_4_ROOT /opt/intel/Compiler/composerxe-2011.4.191)
SET(TBBROOT ${INTEL_12_0_4_ROOT}/tbb)

SET(TBB_INCLUDE_DIRS  ${TBBROOT}/include                                      CACHE FILEPATH "")
SET(TBB_LIBRARY_DIRS  ${TBBROOT}/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21  CACHE FILEPATH "")
