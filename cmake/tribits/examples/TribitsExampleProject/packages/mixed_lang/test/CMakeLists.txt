tribits_add_executable_and_test(
  RayTracerTests
  SOURCES
    tstRay_Tracer.cc
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
