IF (${PACKAGE_NAME}_ENABLE_CUSPARSE)
  # This example depends explicitly and unconditionally on CUSPARSE.
  TRIBITS_ADD_EXAMPLE_DIRECTORIES(hybrid_fem)
ENDIF ()

