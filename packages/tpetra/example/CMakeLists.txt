INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORY(PowerMethod)
ADD_SUBDIRECTORY(HybridPlatform)

IF (Kokkos_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(LocalMatOps)
ENDIF()
