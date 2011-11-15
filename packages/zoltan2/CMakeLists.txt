
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Zoltan2)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
  HAVE_METIS
  "Metis is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PaToH
  HAVE_PATOH
  "PaToH is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ParMETIS
  HAVE_PARMETIS
  "ParMetis is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Scotch
  HAVE_SCOTCH
  "Scotch is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_OVIS
  HAVE_OVIS
  "Ovis is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_HAVE_MALLINFO
  HAVE_MALLINFO
  "malloc.h supports mallinfo()"
  OFF )

#
# The Zoltan library (in zz/zz_util.h) has utilities that
# read /proc/meminfo at runtime.  There is code in Zoltan2
# tests and algorithms that checks this and prints out memory
# currently committed on the compute node.  
# Turning on this option at cmake configure time turns on
# the printing of that info.
#
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LINUX_MEMORY_REPORT
  SHOW_LINUX_MEMINFO
  "turn on runtime display of proc/meminfo data"
  OFF )

#
# The Scotch library keeps track of the maximum amount of memory
# used by Scotch in the memory wrappers in common_memory.c.  But
# there is no function that returns that maximum.  Add this function
# to src/libscotch/common_memory.c:
#
#   size_t SCOTCH_getMemoryMax() { return memorymax; }
#
# and compile scotch with -DCOMMON_MEMORY_TRACE
# The Zoltan2 scotch algorithms will call getMemoryMax and report
# the maximum memory used by Scotch if this cmake option is enabled.
#
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SCOTCH_MEMORY_REPORT
  SHOW_SCOTCH_HIGH_WATER_MARK
  "turn on runtime display maximum memory used by scotch"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_HAVE_TPL64
  HAVE_TPL64
  "64-bit third-party libraries are present"
  OFF )

# 
# Include Fortran support and Fortran tests
#
# IF(Trilinos_ENABLE_Fortran)
# ENDIF()
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
