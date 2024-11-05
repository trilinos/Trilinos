macro(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  assert_defined(${REPOSITORY_NAME}_SOURCE_DIR)
  append_set(CPACK_SOURCE_IGNORE_FILES
    "${${REPOSITORY_NAME}_SOURCE_DIR}/cmake/ctest/"
    )

endmacro()
