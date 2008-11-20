
INCLUDE(Trilinos_Add_Test)

SET(COMMON_SOURCES evilMain.cpp EvilBase.cpp AEvil.cpp BEvil.cpp)

#
# EvilFourFileImplicit
#

TRILINOS_ADD_EXECUTABLE(
  EvilFourFileImplicit
  SOURCES ${COMMON_SOURCES}
  COMM serial mpi  
  )

TRILINOS_ADD_TEST(
  EvilFourFileImplicit
  COMM serial mpi
  )

#
# EvilFourFileExplicit
#

TRILINOS_ADD_EXECUTABLE(
  EvilFourFileExplicit
  SOURCES ${COMMON_SOURCES}
  COMM serial mpi
  )

SET_TARGET_PROPERTIES(
  Teuchos_EvilFourFileExplicit
  PROPERTIES
    COMPILE_DEFINITIONS "DO_EXPLICIT_INSTANTIATION"
  )

TRILINOS_ADD_TEST(
  EvilFourFileExplicit
  COMM serial mpi
  )

# 2008/08/02: rabartl: TODO: The above calls to TRILINOS_ADD_TEST(...) 
# do not find the test executables and I don't know why.
