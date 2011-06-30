
#
# List of extra external repositories that contain extra Trilinos packages.
#
# Extra repositories must be listed here if they are going to be used in
# Continuous (CI) or Nightly testing driven through the
# TrilinosCTestDriverCore.cmake script.
#
# The six fields for each extra external repos are:
#
# NAME: Gives the name of the extra repo (and also the directory name if
# not specified below)
#
# DIR: The relative direcotry location that will be created under the base
# Trilinos base source directory.  This is the name that will be passed to the
# VC tool when cloning/checkingout the repository.
#
# REPOTYPE: The typeof the VC repo (GIT, SVN, or CVS)
#
# REPOURL: This is the URL the extra repository will be cloned (or chekced
# out) from.
#
# PACKSTAT: If the repo has packages, will be "".  If it is NO_PACKAGES, the
# repository does not provide any packges.  A value of NO_PACKAGES is for
# cases where the repo provides code but not add-on Trilinos packages
# (e.g. such as is the case with Dakota used by TriKota).
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

SET( Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY
  preCopyrightTrilinos  ""  GIT  software.sandia.gov:/space/git/preCopyrightTrilinos  ""  Continuous
  TerminalApplication  ""  GIT  software.sandia.gov:/space/git/TerminalApplication  ""   EX 
  )
