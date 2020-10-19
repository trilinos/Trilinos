# Add a single-proc and a multi-proc MueLu test
#
# Usage::
#
#   muelu_add_serial_and_mpi_test( <execRootName>
#     [NAME <testName>]
#     ...
#     NUM_MPI_PROCS <numMpiProcs>
#     ...
#     )
#
# Here, in an MPI build (TPL_ENABLE_MPI=ON), if ``<numMpiProcs> > 1``, then
# two tests will get defined with::
#
#   tribit_add_test( <execRootName>
#     [NAME <testName>]
#     ...
#     NUM_MPI_PROCS 1
#     ...
#     COMM serial mpi
#     )
#
#   tribit_add_test( <execRootName>
#     [NAME <testName>]
#     ...
#     NUM_MPI_PROCS <numMpiProcs>
#     ...
#     COMM mpi
#     )
#
# with the names:
#
# * ``<Package>_<testName>_MPI_1``
# * ``<Package>_<testName>_MPI_<numMpiProcs>``
#
# If ``<numMpiProcs>==1`` in an MPI build, then just the test
# ``<Package>_<testName>_MPI_1`` would get defined.
#
# But in a non-MPI build, the test would get defined as::
#
#   tribit_add_test( <execRootName>
#     [NAME <testName>]
#     ...
#     ...
#     )
#
# and would have the name ``<Package>_<testName>``.  That way, the set of
# tests run in the MPI build are a strick superset of the tests run in a
# non-MPI build.
#
function(muelu_add_serial_and_mpi_test)
  # Get the number of arguments
  list(LENGTH ARGN num_args)
  # Collect the common arguments which is everything excpet for 'NUM_MPI_PROCS
  # <numMpiProcs>' which is captured in the var 'num_mpi_procs')
  set(common_args "")
  set(num_mpi_procs "")
  set(arg_idx 0)
  while(arg_idx LESS ${num_args})
    list(GET ARGN ${arg_idx} arg)
    #print_var(arg)
    if (arg STREQUAL "NUM_MPI_PROCS")
      # Get <numMpiProcs> from the next argument
      math(EXPR arg_idx_p_1 "${arg_idx}+1")
      list(GET ARGN ${arg_idx_p_1} num_mpi_procs)
      math(EXPR arg_idx "${arg_idx}+2")
    else()
      list(APPEND common_args ${arg})
      math(EXPR arg_idx "${arg_idx}+1")
    endif()
  endwhile()
  # Add the test with just 1 proc
  tribits_add_test(${common_args} NUM_MPI_PROCS 1 COMM serial mpi)
  # Add the test with more than 1 proc
  if (${num_mpi_procs} GREATER 1)
    tribits_add_test(${common_args} NUM_MPI_PROCS ${num_mpi_procs} COMM mpi)
  endif()
endfunction()
