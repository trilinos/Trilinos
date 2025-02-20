
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
