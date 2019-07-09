TRIBITS_ADD_EXECUTABLE(
  PinT_TimeStampCoarsening
  SOURCES test_timestamp_coarsening.cpp
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PinT_RestrictProlong
  SOURCES test_restrict_prolong.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 4
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE(
  PinT_HiearchyRebalance
  SOURCES test_hierarchy_coarse_rebalance.cpp
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE(
  PinT_RestrictProlongRebalance
  SOURCES test_restrict_prolong_rebalance.cpp
  ADD_DIR_TO_NAME
  )
