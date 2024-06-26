// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_Utils.h"
#include "Galeri_FiniteElements.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace Galeri;
using namespace Galeri::FiniteElements;

// ==========================================================
// This file solves the scalar problem
//
//   - \mu \nabla u + \sigma u = f    on \Omega
//                           u = g    on \partial \Omega
// 
// where \Omega is a 2D rectangle, divided into triangles.
// The input grid should be generated using similar 
// conventions of file galeri/data/TwoSquares.m:
// - the bc ID of 10 and 20 are for the mortar faces
// - the bc ID of 0 is for the external domain.
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
  if (y < 1)
    return(1.0);
  else
    return(0.5);
}

// Specifies the boundary condition.
int BoundaryType(const int& Patch)
{
  if (Patch == 10 || Patch == 20)
    return(GALERI_DO_NOTHING);
  else
    return(GALERI_DIRICHLET);
}

// Specifies the boundary condition.
double BoundaryValue(const double& x, const double& y, 
                     const double& z, const int& Patch)
{
  if (x == -1.0 || x == 1.0)
    return(0.0);
  else if (Patch == 10 || Patch == 20)
    return(1.0);
  else
    return (0.0);
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

  try {

    // ================================================== //
    // Defines the grid for this problem, a rectangle,    //
    // with the number of nodes along the X-axis (nx) and //
    // Y-axis (ny), the length of the rectangle along the //
    // axes, and the number of processors on each axix.   //
    // ================================================== //

    int nx = 40 * Comm.NumProc();
    int ny = 40;
    int mx = Comm.NumProc();
    int my = 1;

    //TriangleRectangleGrid Grid(Comm, nx, ny, mx, my);
    FileGrid Grid(Comm, "TwoSquares.grid");

    // ======================================================== //
    // Prepares the linear system. This requires the definition //
    // of a quadrature formula compatible with the grid, a      //
    // variational formulation, and a problem object which take //
    // care of filling matrix and right-hand side.              //
    // ======================================================== //
    
    Epetra_CrsMatrix A(Copy, Grid.RowMap(), 0);
    Epetra_Vector    LHS(Grid.RowMap());
    Epetra_Vector    RHS(Grid.RowMap());

    int NumQuadratureNodes = 3;

    GalerkinVariational<TriangleQuadrature>
      Laplace2D(NumQuadratureNodes, Diffusion, Source, Force, 
                BoundaryValue, BoundaryType);

    LinearProblem FiniteElementProblem(Grid, Laplace2D, A, LHS, RHS); 
    FiniteElementProblem.Compute();

    // =================================================== //
    // The solution must be computed here by solving the   //
    // linear system A * LHS = RHS.                        //
    //
    // NOTE: Solve() IS A SIMPLE FUNCTION BASED ON LAPACK, //
    // THEREFORE THE MATRIX IS CONVERTED TO DENSE FORMAT.  //
    // IT WORKS IN SERIAL ONLY.                            //
    // EVEN MEDIUM-SIZED MATRICES MAY REQUIRE A LOT OF     //
    // MEMORY AND CPU-TIME! USERS SHOULD CONSIDER INSTEAD  //
    // AZTECOO, ML, IFPACK OR OTHER SOLVERS.               //
    // =================================================== //

    Solve(&A, &LHS, &RHS);

    // ================== //
    // Output using MEDIT //
    // ================== //
    
    MEDITInterface MEDIT(Comm);
    MEDIT.Write(Grid, "output", LHS);

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
