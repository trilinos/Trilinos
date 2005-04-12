#include "ml_include.h"
#if defined(HAVE_ML_MLAPI)
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "Teuchos_CommandLineProcessor.hpp"
#include "ml_TriangleRectangleGrid.h"
#include "ml_TriangleQuadrature.h"
#include "ml_SUPGVariational.h"
#include "ml_LinearProblem.h"
#include "ml_MEDITInterface.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiLevelNonSymmetricSA.h"
#include "MLAPI_Krylov.h"

// ==========================================================
// Similar to file ml_AdvDiff.cpp, but used MLAPI instead
// of "classical" ML to solve the linear system. The code is
// sligthly more complex than the other examples because 
// the Epetra matrix and vectors have to be wrapped as MLAPI
// object, in order to define an MLAPI preconditioner.
//
// \note For the sake of simplicity, this example can be run
//       with one process only. However, simple changes
//       can make it parallel.
//
// \author Marzio Sala, SNL 9214
//
// \date Last updated on 12-Apr-05.
// ==========================================================

double Diffusion(const double& x, const double& y, const double& z)
{
  return (1.0);
}

double conv = 50000;
double ConvX(const double& x, const double& y, const double& z)
{
  return (conv);
}

double ConvY(const double& x, const double& y, const double& z)
{
  return (-conv);
}

double ConvZ(const double& x, const double& y, const double& z)
{
  return (0.0);
}

double Source(const double& x, const double& y, const double& z)
{
  return (0.0);
}

double Force(const double& x, const double& y, const double& z)
{
  return (0.0);
}

// Specifies the boundary condition.
double BoundaryValue(const double& x, const double& y, 
                     const double& z, const int& Patch)
{
  if ((x == 0.0 && y >= 0.0) || (y == 1.0 && x <= 0.2))
    return(1.0);
  else
    return (0.0);
}

using namespace MLAPI;
using namespace ML_FiniteElements;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // ================================================== //
    // Defines the grid for this problem, a rectangle,    //
    // with the number of nodes along the X-axis (nx) and //
    // Y-axis (ny), the length of the rectangle along the //
    // axes, and the number of processors on each axix.   //
    // ================================================== //

    CommandLineProcessor CLP;

    int nx = 100;
    double DampingFactor = 0.0;
    CLP.setOption("nx", &nx, "number of nodes along the X- and Y-axis");
    CLP.setOption("damp", &DampingFactor, "damping factor");
    CLP.throwExceptions(false);
    CLP.parse(argc,argv);

    if (GetNumProcs() != 1)
    {
      cout << "This example can be run only in serial mode" << endl;
      throw(-1);
    }

    TriangleRectangleGrid Grid(GetEpetra_Comm(), nx, nx, 1, 1);

    // ======================================================== //
    // Prepares the linear system. This requires the definition //
    // of a quadrature formula compatible with the grid, a      //
    // variational formulation, and a problem object which take //
    // care of filling matrix and right-hand side.              //
    // ======================================================== //
    
    Epetra_CrsMatrix A_Epetra(Copy, Grid.RowMap(), 0);
    Epetra_Vector    LHS_Epetra(Grid.RowMap());
    Epetra_Vector    RHS_Epetra(Grid.RowMap());

    int NumQuadratureNodes = 1;

    SUPGVariational<TriangleQuadrature>
      AdvDiff(NumQuadratureNodes, Diffusion, ConvX, 
              ConvY, ConvZ, Source, Force, BoundaryValue);

    LinearProblem FiniteElementProblem(Grid, AdvDiff, A_Epetra, 
                                       LHS_Epetra, RHS_Epetra); 
    FiniteElementProblem.Compute();

    // ======================================================== //
    // Solve the linear system by using MLAPI as preconditioner //
    // ======================================================== //
    
    Teuchos::ParameterList MLList;
    MLList.set("max levels",3);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: damping factor", DampingFactor);
    MLList.set("smoother: type","symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("smoother: damping factor",1.0);
    MLList.set("coarse: max size",32);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type","Amesos-KLU");
    
    Space FineSpace(Grid.RowMap().NumGlobalElements());

    Operator A_MLAPI(FineSpace, FineSpace, &A_Epetra, false);

    MultiVector LHS_MLAPI(FineSpace, LHS_Epetra.Pointers(), 1);
    MultiVector RHS_MLAPI(FineSpace, RHS_Epetra.Pointers(), 1);

    MultiLevelNonSymmetricSA Prec(A_MLAPI, MLList);

    MLList.set("krylov: type", "gmres");
    Krylov(A_MLAPI, LHS_MLAPI, RHS_MLAPI, Prec, MLList);

    // ================== //
    // Output using MEDIT //
    // ================== //
    
    MEDITInterface MEDIT(GetEpetra_Comm());
    MEDIT.Write(Grid, "MLAPI", LHS_Epetra);

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
  puts("--enable-ifpack");
  puts("--enable-amesos");
  puts("--enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif
