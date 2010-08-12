
# This list is just used for unit testing the dependency handling
# CMake code.  The reason that we have a separate list is so that we
# can keep very stable unit tests.
#
# Note: In the real file, the 'DIRS' field must be a single name, not a set of
# directories like it is here.


SET( Trilinos_EXTRAREPOS_DIR_GITREPOURL_CATEGORY
  preCopyrightTrilinos    software.sandia.gov:/git/space/preCopyrightTrilinos   Continuous
  extraTrilinosRepo       software.sandia.gov:/git/space/extraTrilinosRepo      Nightly
  )
