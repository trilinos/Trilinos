INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( PETSC
  REQUIRED_HEADERS petsc.h petscconf.h
  REQUIRED_LIBS_NAMES petscsnes petscksp petscdm petscmat petscvec petsc 
  )
