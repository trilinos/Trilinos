
# -*- cmake -*-

  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E make_directory
      "${Trilinos_BINARY_DIR}/html/common/paramList"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
      "${Trilinos_SOURCE_DIR}/packages/common/paramList"
      "${Trilinos_BINARY_DIR}/html/common/paramList"
    )


