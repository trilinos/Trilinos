tribits_add_test(package1-prg  NOEXEPREFIX  NOEXESUFFIX
  NAME Prg  DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../src"  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Package1 Deps: tpl1" )

tribits_add_advanced_test(Prg-advanced
  TEST_0
    EXEC package1-prg  DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../src"
      NOEXEPREFIX  NOEXESUFFIX
    ARGS "something_extra"
    PASS_REGULAR_EXPRESSION_ALL
      "Package1 Deps: tpl1"
      "something_extra"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  )
