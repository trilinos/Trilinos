
TRIBITS_SUBPACKAGE(Remainder)

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION
  TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION
  "Whether to enable automatic run-time registration of LinearSolverFactory in downstream packages.  All packages that intend to register a LinearSolverFactory at run time, before main() starts, must respect this option.  Its default value is ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION = ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}, because if explicit template instantiation (ETI) is OFF, the automatic registration option increases build time.  The LinearSolver and LinearSolverFactory interfaces remain even if this option is OFF, but in that case, users must do manual registration for all solver packages which they (or solvers) might like to use."
  ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
  )

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_TEST_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()
