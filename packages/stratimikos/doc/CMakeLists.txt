
# -*- cmake -*-

set(DOC_PACKAGE_NAME "stratimikos")

  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E copy
      "${Trilinos_SOURCE_DIR}/packages/${DOC_PACKAGE_NAME}/doc/build_docs"
      "${Trilinos_BINARY_DIR}/packages/${DOC_PACKAGE_NAME}/doc"
    COMMAND "${CMAKE_COMMAND}" -E copy
      "${Trilinos_SOURCE_DIR}/packages/${DOC_PACKAGE_NAME}/doc/build_docs_guts.py"
      "${Trilinos_BINARY_DIR}/packages/${DOC_PACKAGE_NAME}/doc"
    )


#add_subdirectory(parameterList)
