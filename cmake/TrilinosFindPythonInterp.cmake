# Find Python executable which is needed for dependency file building
MACRO(TRILINOS_FIND_PYTHON)
  SET(PythonInterp_FIND_VERSION "2.4")
  ADVANCED_SET(PythonInterp_MUST_BE_FOUND FALSE CACHE BOOL "Require Python to be found or not.") 
  FIND_PACKAGE(PythonInterp)
ENDMACRO()
