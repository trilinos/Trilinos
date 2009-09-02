INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Phalanx)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_EXPLICIT_TEMPLATE_INSTANTIATION
  PHX_ETI
  "Build with explicit template instation enabled."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PHX_DEBUG
  "Build with debug code enabled. Typically, this enables array bounds checking and runtime checking of rank accessors in multidimensional arrays."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_COMPILETIME_ARRAY_CHECK
  PHX_USE_COMPILETIME_ARRAY
  "Enables compile time checking of the rank accessors in the MDArray. This adds safety of checking the rank accessors at compile time.  This can be very inefficient in runtime mode so checking is only done for debug builds if this flag is disabled.  The default is to disable compile time checking.  NOTE: Compile time checking must be OFF for MDFields to be interoperable with Intrepid!"
  OFF )

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
