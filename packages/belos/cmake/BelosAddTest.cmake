function(belos_add_teuchos_and_kokkos_test)
  # Get the number of arguments
  list(LENGTH ARGN num_args)
  # Collect the common arguments which is everything except for 'NAME'
  set(common_args "")
  set(arg_idx 0)
  while(arg_idx LESS ${num_args})
    list(GET ARGN ${arg_idx} arg)
    #print_var(arg)
    if (arg STREQUAL "NAME")
      # Get <NAME> from the next argument
      math(EXPR arg_idx_p_1 "${arg_idx}+1")
      list(GET ARGN ${arg_idx_p_1} name)
      math(EXPR arg_idx "${arg_idx}+2")
    else()
      list(APPEND common_args ${arg})
      math(EXPR arg_idx "${arg_idx}+1")
    endif()
  endwhile()
  # Add the test with the Teuchos environment variable
  tribits_add_test(${common_args} NAME ${name}_Teuchos ENVIRONMENT BELOS_DENSE_MATRIX_ABSTRACTION=Teuchos)
  # Add the test with the Kokkos environment variable
  tribits_add_test(${common_args} NAME ${name}_Kokkos ENVIRONMENT BELOS_DENSE_MATRIX_ABSTRACTION=Kokkos)
endfunction()


function(belos_add_executable_and_teuchos_and_kokkos_tests exe)
  # Get the number of arguments
  list(LENGTH ARGN num_args)
  # Collect the common arguments which is everything except for 'SOURCES'
  set(common_args "")
  set(arg_idx 0)
  while(arg_idx LESS ${num_args})
    list(GET ARGN ${arg_idx} arg)
    #print_var(arg)
    if (arg STREQUAL "SOURCES")
      # Get <NAME> from the next argument
      math(EXPR arg_idx_p_1 "${arg_idx}+1")
      list(GET ARGN ${arg_idx_p_1} sources)
      math(EXPR arg_idx "${arg_idx}+2")
    else()
      list(APPEND common_args ${arg})
      math(EXPR arg_idx "${arg_idx}+1")
    endif()
  endwhile()
  # Add the executable
  tribits_add_executable(${exe} SOURCES ${sources})
  # Add the tests
  belos_add_teuchos_and_kokkos_test(${exe} NAME ${exe} ${common_args})
endfunction()
