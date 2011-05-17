
#
# List of extra external repositories that contain extra Trilinos packages.
#
# Extra repositories must be listed here if they are going to be used in
# Continuous (CI) or Nightly testing driven through the
# TrilinosCTestDriverCore.cmake script.
#
# The three fields for each extra external repos are:
#
# DIR: Gives the name of the extra repo that will be created under the base
# Trilinos base source directory.  This needs to be a single name.  This is
# the name that will be passed to git when cloning the repository.
#
# GITREPOURL: This is the git repo URL where the extra repository will be
# cloned from.
#
# CATEGORY: This is category of tests where extra repo will be pulled in for.
# Valid categories include:
#
#   Continuous: Continuous integration testing run throughout the development day
#
#   Nightly: Nightly testing run at the end of every development day
#
#   EX: Experimental, not run implicitly for any type of testing, including
#     for 'Experimental' builds.
#
# The reason that you would want to list 'EX' repositories, is that it leaves
# the option to explicitly add the repositories in specialized nightly
# testing.
#
# NOTE: The extra repositories must be listed in assending order of package
# dependencies, just like Trilinos packages in the file
# TrilinosPackages.cmake.  This also means that there can not be any cicular
# dependencies between groups of packages within an extra repository.
#
# NOTE: The packages in a downstream extra repo can depend on packages in an
# upstream repo.
#

SET( Trilinos_EXTRAREPOS_DIR_GITREPOURL_CATEGORY
  LIMEExt  software.sandia.gov:/space/git/LIMEExt  Continuous
  PSSDriversExt  nstdsrv.ornl.gov/gitroot/casl_pssdrivers  Continuous
  )
