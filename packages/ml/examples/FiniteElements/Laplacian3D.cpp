#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "AztecOO.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_AbstractGrid.h"
#include "ml_HexCubeGrid.h"
#include "ml_HexQuadrature.h"
#include "ml_GalerkinVariational.h"
#include "ml_LinearProblem.h"
#include "ml_MEDITInterface.h"

// ==========================================================
// This file solves the scalar elliptic problem
//
//   - \mu \nabla u + \sigma u = f    on \Omega
//                           u = g    on \partial \Omega
// 
// where \Omega is a 3D cube, divided into hexahedra.
// `f' is specified by function `Force()', the Dirichlet boundary condition
// by function `BoundaryValue()', and the value of \mu and
// \sigma can be changed in the functions Diffusion() and
// Source(). The code solves the corresponding linear system 
// using AztecOO with ML preconditioner, and writes the 
// solution into a MEDIT-compatible format. Then, it computes 
// the L2 and H1 norms of the solution and the error.
//
// \note This example requires the number of processors to
//       be the cube of an integer number.
//
// \author Marzio Sala, SNL 9214
//
// \date Last updated on 31-Mar-05.
// ==========================================================

double Diffusion(const double& x, const double& y, const double& z)
{
  return (1.0);
}

double Source(const double& x, const double& y, const double& z)
{
  return (0.0);
}

double Force(const double& x, const double& y, const double& z)
{
  return (-6.0);
}

// Specifies the boundary condition.
double BoundaryValue(const double& x, const double& y, 
                     const double& z, const int& Patch)
{
  return(x * x + y * y + z * z);
}

// Returns the value of the exact solution and its first
// derivatives with respect to x, y and z.
int ExactSolution(double x, double y, double z, double* res)
{
  res[0] = x * x + y * y + z * z;
  res[1] = 2 * x;
  res[2] = 2 * y;
  res[3] = 2 * z;

  return(0);
}

using namespace ML_Epetra;
using namespace ML_FiniteElements;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  try {

    // ================================================== //
    // Defines the grid for this problem, a rectangle,    //
    // with the number of nodes along the X-axis (nx) and //
    // Y-axis (ny), the length of the rectangle along the //
    // axes, and the number of processors on each axix.   //
    // ================================================== //

    CommandLineProcessor CLP;

    int nx = 10;
    int ny = 10;
    int nz = 10;
    int mx = (int)pow((double)Comm.NumProc(), 0.3334);
    int my = mx;
    int mz = mx;
    CLP.setOption("nx", &nx, "number of nodes along the X-axis");
    CLP.setOption("ny", &ny, "number of nodes along the Y-axis");
    CLP.setOption("nz", &nz, "number of nodes along the Z-axis");
    CLP.setOption("mx", &mx, "number of subdomains along the X-axis");
    CLP.setOption("my", &my, "number of subdomains along the Y-axis");
    CLP.setOption("mz", &mz, "number of subdomains along the Z-axis");

    CLP.throwExceptions(false);
    CLP.parse(argc,argv);

    if (mx * my * mz != Comm.NumProc()) {
      if (Comm.MyPID() == 0) 
      {
        cout << "Number of subdomains, " << mx << " x " << my << " x " << mz << endl;
        cout << "does not equal to the total number of processes (" << Comm.NumProc() << endl;
        cout << "Please re-run with --help option for details." << endl;
      }
      throw(-1);
    }

    HexCubeGrid Grid(Comm, nx, ny, nz, mx, my, mz);

    // ======================================================== //
    // Prepares the linear system. This requires the definition //
    // of a quadrature formula compatible with the grid, a      //
    // variational formulation, and a problem object which take //
    // care of filling matrix and right-hand side.              //
    // ======================================================== //
    
    Epetra_CrsMatrix A(Copy, Grid.RowMap(), 0);
    Epetra_Vector    LHS(Grid.RowMap());
    Epetra_Vector    RHS(Grid.RowMap());

    int NumQuadratureNodes = 1;

    GalerkinVariational<HexQuadrature>
      Laplacian(NumQuadratureNodes, Diffusion, Source, Force, 
                BoundaryValue);

    LinearProblem FiniteElementProblem(Grid, Laplacian, A, LHS, RHS); 
    FiniteElementProblem.Compute();

    // ======================================================= //
    // Solve the linear system by using ML as a preconditioner //
    // and AztecOO's conjugate gradient as a solver.           //
    // ======================================================= //
    
    Teuchos::ParameterList MLList;
    MLList.set("max levels",3);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: damping factor", 1.333);
    MLList.set("smoother: type","symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("smoother: damping factor",1.0);
    MLList.set("coarse: max size",32);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type","Amesos-KLU");
    
    MultiLevelPreconditioner MLPrec(A, MLList);

    Epetra_LinearProblem Problem(&A, &LHS, &RHS);
    AztecOO Solver(Problem);
    Solver.SetPrecOperator(&MLPrec);

    Solver.SetAztecOption(AZ_solver, AZ_cg);
    Solver.SetAztecOption(AZ_output, 32);
    Solver.Iterate(1550, 1e-10);

    FiniteElementProblem.ComputeNorms(LHS, ExactSolution);
    
    // ================== //
    // Output using MEDIT //
    // ================== //
    
    MEDITInterface MEDIT(Comm);
    MEDIT.Write(Grid, "Laplacian3D", LHS);

  }
  catch (int e) {
    cerr << "Caught exception, value = " << e << endl;
  }
  catch (...) {
    cerr << "Caught generic exception" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-amesos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif
