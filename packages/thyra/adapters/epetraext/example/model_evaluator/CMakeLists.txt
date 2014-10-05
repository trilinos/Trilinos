INCLUDE(TribitsAddExecutableAndTest)


INCLUDE_DIRECTORIES(
  ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/2dsim
  ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/4dopt
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ForwardSolveEpetraModelEval2DSim
  SOURCES
    ForwardSolveEpetraModelEval2DSimMain.cpp
    ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/2dsim/EpetraModelEval2DSim.cpp
    ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/4dopt/EpetraModelEval4DOpt.cpp
  COMM serial mpi
  )
