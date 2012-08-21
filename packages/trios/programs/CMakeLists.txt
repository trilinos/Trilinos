INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(programs)

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

############# Special exception for rpcgen extra commas in enum ##########
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

  TRIBITS_ADD_EXECUTABLE(
    kill-service
    SOURCES nssi-kill.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )


  TRIBITS_ADD_EXECUTABLE(
    ping-service
    SOURCES nssi-ping.cpp
    DEPLIBS ${DEPLIBS}
    NOEXEPREFIX
  )

TRIBITS_SUBPACKAGE_POSTPROCESS()
