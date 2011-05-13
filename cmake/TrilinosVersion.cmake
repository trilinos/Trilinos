#
# Single file that needs to be changed on a release branch
# or on the development branch in order to configure Trilinos
# for release mode and set the version.
#

SET(Trilinos_VERSION "10.7")
SET(Trilinos_MAJOR_VERSION "10")
SET(Trilinos_MAJOR_MINOR_VERSION "100700")
SET(Trilinos_VERSION_STRING "10.7 (Dev)")
SET(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT ON) # Change to 'OFF' for a release

# Used by testing scripts and should not be used elsewhere
SET(Trilinos_REPOSITORY_BRANCH "master" CACHE INTERNAL "")
SET(Trilinos_TESTING_TRACK "" CACHE INTERNAL "")

SET(TRILINOS_MAJOR_VERSION TRUE)
SET(TRILINOS_MAJOR_MINOR_VERSION TRUE)
SET(TRILINOS_VERSION_STRING TRUE)
