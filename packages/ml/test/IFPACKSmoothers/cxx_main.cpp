#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

using namespace Trilinos_Util;
using namespace Teuchos;
using namespace ML_Epetra;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // ===================== //
  // create linear problem //
  // ===================== //
  
  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", 27000);
  
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_MultiVector& LHS = *(Gallery.GetStartingSolution());
  Epetra_MultiVector& RHS = *(Gallery.GetRHS());

  AztecOO solver(*Problem);

  bool TestPassed = true;

  vector<string> TestList;
  TestList.push_back("Jacobi");
  TestList.push_back("Gauss-Seidel");
  TestList.push_back("symmetric Gauss-Seidel");

  double Damping = 0.67;

  for (int sweeps = 1 ; sweeps < 5 ; sweeps += 2) {

    for (unsigned int i = 0 ; i < TestList.size() ; ++i) {

      PrintLine();
      cout << "### Testing " << TestList[i] << endl;
      cout << "### sweeps = " << sweeps << endl;
      cout << "### damping = " << Damping << endl;
      PrintLine();

      // ========================= //
      // build ML with ML smoother //
      // ========================= //

      ParameterList MLList;
      SetDefaults("SA",MLList);
      MLList.set("smoother: type", TestList[i]);
      MLList.set("smoother: pre or post", "both");
             MLList.set("smoother: type", TestList[i]);
      MLList.set("smoother: sweeps", sweeps);
      MLList.set("smoother: damping factor", Damping);

      MLList.set("output", 0);

      Epetra_Time Time(Comm);

      MultiLevelPreconditioner* MLPrec = 
        new MultiLevelPreconditioner(*A, MLList);

      // tell AztecOO to use this preconditioner, then solve
      solver.SetPrecOperator(MLPrec);

      solver.SetAztecOption(AZ_solver, AZ_gmres);
      solver.SetAztecOption(AZ_output, 32);

      LHS.PutScalar(0.0);
      RHS.PutScalar(1.0);

      solver.Iterate(1550, 1e-12);

      delete MLPrec;
      MLPrec = 0;

      int MLIters = solver.NumIters();
      double MLTime = Time.ElapsedTime();

      // ============================= //
      // build ML with IFPACK smoother //
      // ============================= //

      // reset the options, and stick IFPACK parameters list as 
      // needed by ML
      SetDefaults("SA",MLList);
      MLList.set("smoother: pre or post", "both");
      MLList.set("smoother: type", "IFPACK");
      MLList.set("smoother: ifpack type", "point relaxation (no AS)");
      ParameterList& IFPACKList = MLList.sublist("smoother: ifpack list");;
      IFPACKList.set("relaxation: type", TestList[i]);
      IFPACKList.set("relaxation: sweeps", sweeps);
      IFPACKList.set("relaxation: damping factor", Damping);
      IFPACKList.set("relaxation: zero starting solution", false);
      IFPACKList.set("partitioner: type", "metis");
      IFPACKList.set("partitioner: local parts", -1);

      MLList.set("output", 0);
      Time.ResetStartTime();

      MLPrec = new MultiLevelPreconditioner(*A, MLList);

      // tell AztecOO to use this preconditioner, then solve
      solver.SetPrecOperator(MLPrec);

      solver.SetAztecOption(AZ_solver, AZ_gmres);
      solver.SetAztecOption(AZ_output, 32);

      LHS.PutScalar(0.0);
      RHS.PutScalar(1.0);

      solver.Iterate(1550, 1e-12);

      delete MLPrec;

      int IFPACKIters = solver.NumIters();
      double IFPACKTime = Time.ElapsedTime();
      Time.ResetStartTime();

      // check whether we get the same results or not,
      // and also compare the CPU time for both
      if (MLIters != IFPACKIters) {
        cout << "TEST FAILED: ML converged in " << MLIters << ", while" << endl;
        cout << "IFPACK converged in " << IFPACKIters << " iterations." << endl;
        TestPassed = false;
      }

      if (Comm.MyPID() == 0) {
        cout << "ML time     = " << MLTime << " (s)" << endl;
        cout << "IFPACK time = " << IFPACKTime << " (s)" << endl;
      }
    }
  }

  if (!TestPassed && Comm.MyPID() == 0) {
    cerr << "Test `IFPACKSmoothers.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0) 
    cout << "Test `IFPACKSmoothers.exe' passed!" << endl;
  return (EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
    
  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils --enable-ifpack");

#ifdef HAVEML_MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
