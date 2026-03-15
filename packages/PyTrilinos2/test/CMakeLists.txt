
ENABLE_TESTING()

INCLUDE(PyTrilinos2MakeTest)

PyTrilinos2_MAKE_MPI_TEST(CG)
PyTrilinos2_MAKE_MPI_TEST(parameterList)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Stratimiko_cp
    SOURCE_FILES Stratimikos.py)

TRIBITS_ADD_TEST(
  ${Python3_EXECUTABLE}
  NOEXEPREFIX
  NOEXESUFFIX
  NAME Stratimikos
  POSTFIX_AND_ARGS_0 "LU"                 Stratimikos.py --problemSize=100   --solver=LU
  POSTFIX_AND_ARGS_1 "CG"                 Stratimikos.py --problemSize=100   --solver=CG       --prec=None
  POSTFIX_AND_ARGS_2 "CG_Jacobi"          Stratimikos.py --problemSize=100   --solver=CG       --prec=Jacobi
  POSTFIX_AND_ARGS_3 "BiCGStab_Chebyshev" Stratimikos.py --problemSize=100   --solver=BiCGStab --prec=Chebyshev
  POSTFIX_AND_ARGS_4 "GMRES_ILU"          Stratimikos.py --problemSize=1000  --solver=GMRES    --prec=ILU
  POSTFIX_AND_ARGS_5 "CG_multigrid"       Stratimikos.py --problemSize=10000 --solver=CG       --prec=multigrid
  PASS_REGULAR_EXPRESSION "OK"
  ENVIRONMENT "PYTHONPATH=${PyTrilinos2_PYTHONPATH}:$ENV{PYTHONPATH}"
)

find_package(pybind11 REQUIRED)
find_package(Torch REQUIRED)
MESSAGE("TORCH_CXX_FLAGS: ${TORCH_CXX_FLAGS}")
MESSAGE("TORCH_LIBRARIES: ${TORCH_LIBRARIES}")

find_library(TORCH_PYTHON_LIBRARY torch_python PATH "${TORCH_INSTALL_PREFIX}/lib")
MESSAGE("TORCH_PYTHON_LIBRARY: ${TORCH_PYTHON_LIBRARY}")


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RTC
  SOURCES rtc.cpp
  ENVIRONMENT "PYTHONPATH=${PyTrilinos2_PYTHONPATH}:$ENV{PYTHONPATH}"
  NUM_MPI_PROCS 1
)

target_link_libraries(PyTrilinos2_RTC PRIVATE pybind11::embed)
target_link_libraries(PyTrilinos2_RTC PRIVATE ${TORCH_LIBRARIES} ${TORCH_PYTHON_LIBRARY})
