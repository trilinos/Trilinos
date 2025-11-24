#
# Single file that needs to be changed on a release branch
# or on the development branch in order to configure Trilinos
# for release mode and set the version.
#

SET(Trilinos_VERSION 16.2.0)
SET(Trilinos_MAJOR_VERSION 16)
SET(Trilinos_MAJOR_MINOR_VERSION 160200)
SET(Trilinos_VERSION_STRING "16.2.0")
SET(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT OFF) # Change to 'OFF' for a release

# Used by testing scripts and should not be used elsewhere
SET(Trilinos_REPOSITORY_BRANCH "trilinos-release-16-2-branch" CACHE INTERNAL "")
SET(Trilinos_EXTRAREPOS_BRANCH "master" CACHE INTERNAL "")
SET(Trilinos_TESTING_TRACK "" CACHE INTERNAL "")

# NOTE: Above, the extra repos for Trilinos don't have a 'develop' branch yet
# so you have to run this with 'master'.  But on a release branch, these
# should be change to all the same branch.
