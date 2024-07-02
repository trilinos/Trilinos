
#This find is largely not useful. The blacs libraries by default do not use a standard
#naming convention. They are called "blacs.a" instead of "libblacs.a" The TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES
#currently cannot find them like that so you must specify them with 
#TPL_SCALAPACK_LIBRARIES=<path>/libscalapack.a;<path>/blacs.a;<path>blacsF77.a;<path>blacs.a
#hopefully this will eventually be fixed, either in blacs or a workaround done in cmake.

#The blacs libraries apparently have a circular dependency between blacs and blacsF77
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SCALAPACK
  REQUIRED_LIBS_NAMES scalapack blacs blacsF77 blacs)
