
# We need to inject the Trilinos/cmake directory to find
# TrilinosCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Trilinos_SOURCE_DIR}/cmake")


MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  #MESSAGE("TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS got called!")

  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )
    
  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting Trilinos_ENABLE_ForTrilinos=OFF"
      " because Trilinos_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(Trilinos_ENABLE_ForTrilinos OFF)
  ENDIF()
    
  # ToDo: What is this and why is it needed?
  SET(TRILINOS_BUILD_SHARED_LIBS "@BUILD_SHARED_LIBS@")

ENDMACRO()
