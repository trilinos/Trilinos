function(belos_add_teuchos_and_kokkos_test)

  # determine how many run lines got passed via ARGS
  cmake_parse_arguments(
    PARSE_ARGV 1 # One named argument
    #prefix
    PARSE
    # options
    "NOEXEPREFIX;NOEXESUFFIX;STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;RUN_SERIAL"
    #one_value_keywords
    "DISABLED;LIST_SEPARATOR"
    #multi_value_keywords
    "DIRECTORY;KEYWORDS;COMM;NUM_MPI_PROCS;NUM_TOTAL_CORES_USED;ARGS;${POSTFIX_AND_ARGS_LIST};NAME;NAME_POSTFIX;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;EXCLUDE_IF_NOT_TRUE;PASS_REGULAR_EXPRESSION;FAIL_REGULAR_EXPRESSION;TIMEOUT;ENVIRONMENT;ADDED_TESTS_NAMES_OUT"
  )
  list(LENGTH PARSE_ARGS num_arg_args)

  # Get the number of arguments
  list(LENGTH ARGN num_args)
  # Collect the common arguments which is everything except for 'NAME' and 'ARGS'
  set(common_args "")
  set(common_arg_args "")
  set(arg_idx 0)
  while(arg_idx LESS ${num_args})
    list(GET ARGN ${arg_idx} arg)

    if (arg STREQUAL "NAME")
      # Get <NAME> from the next argument
      math(EXPR arg_idx_p_1 "${arg_idx}+1")
      list(GET ARGN ${arg_idx_p_1} name)
      math(EXPR arg_idx "${arg_idx}+2")
    elseif(arg STREQUAL "ARGS")
      # Get all ARGS
      set(arg_arg_idx 0)
      while(arg_arg_idx LESS ${num_arg_args})
        math(EXPR arg_idx_p_1 "${arg_idx}+1+${arg_arg_idx}")
        list(GET ARGN ${arg_idx_p_1} arg_arg)
        list(APPEND common_arg_args ${arg_arg})
        math(EXPR arg_arg_idx "${arg_arg_idx}+1")
      endwhile()
      math(EXPR arg_idx "${arg_idx}+1+${num_arg_args}")
    else()
      list(APPEND common_args ${arg})
      math(EXPR arg_idx "${arg_idx}+1")
    endif()
  endwhile()

  # Option values to be passed as commandline args
  set(denseMatrices "")
  list(APPEND denseMatrices "Teuchos")
  list(APPEND denseMatrices "Kokkos")

  # Add tests with augmented ARGS and NAME
  foreach(denseMatrix IN LISTS denseMatrices)
    set(processed_arg_args "")
    foreach(arg IN LISTS common_arg_args)
      list(APPEND processed_arg_args "${arg} --denseMatrix=${denseMatrix}")
    endforeach()

    tribits_add_test(${common_args} NAME ${name}_${denseMatrix} ARGS ${processed_arg_args})
  endforeach()


  # Add the test with the Teuchos environment variable
  # tribits_add_test(${common_args} NAME ${name}_Teuchos ENVIRONMENT BELOS_DENSE_MATRIX_ABSTRACTION=Teuchos ARGS ${common_arg_args})
  # Add the test with the Kokkos environment variable
  # tribits_add_test(${common_args} NAME ${name}_Kokkos ENVIRONMENT BELOS_DENSE_MATRIX_ABSTRACTION=Kokkos ARGS ${common_arg_args})
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
