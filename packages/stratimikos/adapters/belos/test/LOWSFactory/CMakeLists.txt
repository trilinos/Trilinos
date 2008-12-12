INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE(
  test_single_belos_thyra_solver_driver
  SOURCES
    test_single_belos_thyra_solver_driver.cpp
    test_single_belos_thyra_solver.cpp
  COMM serial mpi
  )


TRILINOS_ADD_TEST(
  test_single_belos_thyra_solver_driver
  ARGS
    "--max-iters=4 --no-use-preconditioner --show-all-tests  --matrix-file=FourByFour.mtx"
    "--show-all-tests --gmres-krylov-length=190 --max-solution-error=1e-2 --max-resid=1e-6 --matrix-file=nos5.mtx"
    "--no-use-preconditioner --show-all-tests --num-rhs=1 --block-size=1 --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6 --matrix-file=nos1.mtx"
    "--no-use-preconditioner --show-all-tests --num-rhs=8 --block-size=8 --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6 --matrix-file=nos1.mtx"
    "--show-all-tests --num-rhs=1 --block-size=1 --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6 --matrix-file=nos1.mtx"
    "--show-all-tests --num-rhs=8 --block-size=8 --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6 --matrix-file=nos1.mtx"
  COMM serial
  )


COPY_FILES_TO_BINARY_DIR(CopyBelosTestFiles
  DEST_FILES
    bcsstk01.mtx                              
    bcsstk02.mtx                              
    bcsstk04.mtx                              
    nos1.mtx                                      
    nos5.mtx                                      
    FourByFour.mtx                          
    Diagonal.mtx                              
    MissingADiagonal.mtx                              
    KheadK.mtx                                  
    KheadSorted.mtx
    Khead.triS  
  SOURCE_PREFIX
    "In_"
  SOURCE_DIR
    "${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic"
  EXEDEPS
    test_single_belos_thyra_solver_driver
  )
