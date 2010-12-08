# Required modules

INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)

# Example definitions

ADD_SUBDIRECTORY(C)

IF (Trilinos_ENABLE_CXX)
ADD_SUBDIRECTORY(CPP)
ENDIF()

#IF (Trilinos_ENABLE_Fortran AND Zoltan_ENABLE_F90INTERFACE)
#ADD_SUBDIRECTORY(Fortran)
#ENDIF()
