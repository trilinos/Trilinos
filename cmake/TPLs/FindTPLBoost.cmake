INCLUDE(Assert_Defined)
INCLUDE(Global_Null_Set)
INCLUDE(Advanced_Set)


IF (TPL_Boost_INCLUDE_DIRS)

  # The Boost library has already been found or the user has specified
  # these manually.  In this case, just make sure that everything has been
  # specified correctly

  # Make sure the other variables have been defined
  ASSERT_DEFINED(TPL_Boost_LIBRARY_DIRS)
  ASSERT_DEFINED(TPL_Boost_LIBRARIES)

  # Verify that indeed we have found Boost!

  # ToDo: Implement!

ELSE()

  # Otherwise, we need to look for the Boost headers and libraries

  # This is the CMake built-in FindBoost module
  FIND_PACKAGE(Boost COMPONENTS iostreams)

  IF (Boost_FOUND)

     ASSERT_DEFINED(Boost_INCLUDE_DIRS)
     ADVANCED_SET(TPL_Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} CACHE STRING
       "Boost include directories" )

     ASSERT_DEFINED(Boost_LIBRARY_DIRS)
     ADVANCED_SET(TPL_Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIRS} CACHE STRING
       "Boost library directories" )

     ADVANCED_SET(TPL_Boost_LIBRARIES ${Boost_LIBRARIES} CACHE STRING
       "Boost library directories" )
     # ToDo: Add whatever boost libraries you want above!  Or, the user
     # can add them!

  ELSE()

    MESSAGE(FATAL_ERROR "Error, could not find the boost Library!")

  ENDIF()

ENDIF()
