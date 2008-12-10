
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

TRILINOS_PACKAGE(Stratimikos)

ADD_SUBDIRECTORY(adapters)

ADD_SUBDIRECTORY(src)

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

TRILINOS_PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

TRILINOS_PACKAGE_POSTPROCESS()

#
# ToDos:
#
# (*) 2008/12/10: rabartl: Add support for Belos once Belos has been
# CMakeifyied
#
# (*) 2008/12/10: rabartl: Add general support for generating sample
# output used by doxygen documentation.
#
# (*) 2008/12/10: rabartl: Change TRILINOS_PACKAGE_ADD_LIBRARY(...) to
# accept DEPPACKAGES.  Then, use this with the adatper and main
# libraries in Stratimikos to better define the library dependency
# structure.  Right now, for example, the thyra/amesos adatpers depend
# on ML, AztecOO, etc when they should not.
#
