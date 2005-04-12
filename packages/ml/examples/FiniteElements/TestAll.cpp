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
#include "ml_TriangleRectangleGrid.h"
#include "ml_TriangleQuadrature.h"
#include "ml_QuadRectangleGrid.h"
#include "ml_QuadQuadrature.h"
#include "ml_TetCubeGrid.h"
#include "ml_TetQuadrature.h"
#include "ml_HexCubeGrid.h"
#include "ml_HexQuadrature.h"
#include "ml_GalerkinVariational.h"
#include "ml_LinearProblem.h"

// ==========================================================
// \brief Simple test for all grid and quadrature classes.
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
  return (0.0);
}

double BoundaryValue(const double& x, const double& y, 
                     const double& z, const int& Patch)
{
  return(x + y);
}

int ExactSolution(double x, double y, double z, double* res)
{
  res[0] = x + y;
  res[1] = x;
  res[2] = y;
  res[3] = 0.0;

  return(0);
}

using namespace ML_Epetra;
using namespace ML_FiniteElements;

// ====================================================================== 
// Test hexahedral elements.
// ====================================================================== 

template<class T>
bool Test(const AbstractGrid& Grid)
{
  int NumQuadratureNodes = 1;
  GalerkinVariational<T> Laplacian(NumQuadratureNodes, Diffusion, Source, 
                                   Force, BoundaryValue);

  Epetra_CrsMatrix A(Copy, Grid.RowMap(), 0);
  Epetra_Vector    LHS(Grid.RowMap());
  Epetra_Vector    RHS(Grid.RowMap());

  LinearProblem FiniteElementProblem(Grid, Laplacian, A, LHS, RHS); 
  FiniteElementProblem.Compute();

  Teuchos::ParameterList MLList;
  MLList.set("output", 0);

  MultiLevelPreconditioner MLPrec(A, MLList);

  Epetra_LinearProblem Problem(&A, &LHS, &RHS);
  AztecOO Solver(Problem);
  Solver.SetPrecOperator(&MLPrec);

  Solver.SetAztecOption(AZ_solver, AZ_cg);
  Solver.SetAztecOption(AZ_output, AZ_none);
  Solver.Iterate(1550, 1e-10);

  double SolutionNorm[3];
  double ExactNorm[3];
  double DiffNorm[3];

  FiniteElementProblem.ComputeNorms(LHS, ExactSolution, false, SolutionNorm, 
                                    ExactNorm, DiffNorm);
  
  if (DiffNorm[0] > 1.0e-4)
    return(false);
  else
    return(true);
}

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

  int n = 10;

  if (Comm.NumProc() != 1)
  {
    cerr << "This test is only serial" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  if (!Test<TriangleQuadrature>(TriangleRectangleGrid(Comm, n, n, 1, 1))) 
    cerr << "Test with triangles failed!" << endl;
  else
    cerr << "Test with triangles passed!" << endl;

  if (!Test<QuadQuadrature>(QuadRectangleGrid(Comm, n, n, 1, 1))) 
    cerr << "Test with quads failed!" << endl;
  else
    cerr << "Test with quads passed!" << endl;

  if (!Test<TetQuadrature>(TetCubeGrid(Comm, n, n, n, 1, 1, 1))) 
    cerr << "Test with tetrahedra failed!" << endl;
  else
    cerr << "Test with tetrahedra passed!" << endl;

  if (!Test<HexQuadrature>(HexCubeGrid(Comm, n, n, n, 1, 1, 1))) 
    cerr << "Test with hexahedra failed!" << endl;
  else
    cerr << "Test with hexahedra passed!" << endl;

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
