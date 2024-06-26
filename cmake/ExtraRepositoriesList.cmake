#
# List of extra repositories that contain extra Trilinos packages.
#
# See documentation for the function
# TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES() in the TriBITS Guide and
# Reference:
#
#  https://tribits.org/doc/TribitsDevelopersGuide.html#tribits-project-define-extra-repositories
#

IF (NOT "$ENV{Trilinos_REPOS_URL_BASE}" STREQUAL "")
  SET(Trilinos_REPOS_URL_BASE  $ENV{Trilinos_REPOS_URL_BASE}
    CACHE STRING "Set from env var Trilinos_REPOS_URL_BASE" FORCE)
ELSE()
  SET(Trilinos_REPOS_URL_BASE_DEFAULT git@github.com:trilinos/)
  SET(Trilinos_REPOS_URL_BASE  ${Trilinos_REPOS_URL_BASE_DEFAULT}
    CACHE STRING "Base URL to Trilinos repos <url-base><repo-name>")
ENDIF()
MARK_AS_ADVANCED(Trilinos_REPOS_URL_BASE)

TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  MOOCHO_repo  packages/moocho  GIT  ${Trilinos_REPOS_URL_BASE}moocho  NOPACKAGES  Nightly
  CTrilinos_repo  packages/CTrilinos  GIT  ${Trilinos_REPOS_URL_BASE}CTrilinos  NOPACKAGES  Nightly
  Optika_repo  packages/optika  GIT  ${Trilinos_REPOS_URL_BASE}optika  NOPACKAGES  Nightly
  Mesquite_repo  packages/mesquite  GIT  ${Trilinos_REPOS_URL_BASE}mesquite  NOPACKAGES  Nightly
  Avatar_repo packages/avatar  GIT  https://github.com/sandialabs/avatar.git  NOPACKAGES Nightly
  preCopyrightTrilinos  ""  GIT  software.sandia.gov:/space/git/preCopyrightTrilinos  ""  Continuous
  TerminalApplication  ""  GIT  software.sandia.gov:/space/git/TerminalApplication  ""   EX
  )
