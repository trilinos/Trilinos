INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt) 

    TRIBITS_ADD_EXECUTABLE(
      CrashMe1
      SOURCES crashme1.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )
     
    #TRIBITS_ADD_EXECUTABLE(
    #  CrashMe
    #  SOURCES visualizeAggregates2.cpp
    #  COMM serial mpi
    #  DEPLIBS muelu xpetra
    #  )
 
    # TRIBITS_ADD_EXECUTABLE(
    #  CrashMe2
    #  SOURCES visualizeAggregates3.cpp
    #  COMM serial mpi
    #  DEPLIBS muelu xpetra
    #  )
    
    TRIBITS_COPY_FILES_TO_BINARY_DIR(CrashMe_example_cp
      SOURCE_FILES Condif2Mat.mat Condif2Rhs.mat
      )


ENDIF()
