INCLUDE(TPLDeclareLibraries)

# Thrust TPL requires CUDA
IF (NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "\nThrust TPL requires that CUDA support is enabled. Please set \n  TPL_ENABLE_CUDA=ON\n\n")
ELSE()
  TPL_DECLARE_LIBRARIES( Thrust
      REQUIRED_HEADERS thrust/for_each.h
  )
ENDIF()

