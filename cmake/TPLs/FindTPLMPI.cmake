INCLUDE(TPLDeclareLibraries)

# Either the MPI compiler wrappers take care of these or the user has to set
# the explicitly using basic compile flags and Trilinos_EXTRA_LINK_FLAGS.
GLOBAL_SET(TPL_MPI_INCLUDE_DIRS)
GLOBAL_SET(TPL_MPI_LIBRARIES)
GLOBAL_SET(TPL_MPI_LIBRARY_DIRS)
