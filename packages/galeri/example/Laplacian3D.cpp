// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#include "Galeri_Utils.h"
#include "Galeri_FiniteElements.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

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
// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 15-Sep-05
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
                     const double& z, const int& PatchID)
{
  return(x * x + y * y + z * z);
}

// Specifies the boundary condition.
int BoundaryType(const int& PatchID)
{
  return(Galeri::FiniteElements::GALERI_DIRICHLET);
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

using namespace Galeri;
using namespace Galeri::FiniteElements;

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

    // ============================================================ //
    // Prepares the computational domain. For simplicity,           //
    // the computation domain has (nx * NumProcs, ny, nz) elements, //
    // and it is partitioned into (NumProcs, 1, 1) subdomains.      //
    // If you want to change the grid element, remember to modify   //
    // the quadrature in GalerkinVariational<T>. Now T is set to    //
    // HexQuadrature.                                               //
    // ============================================================ //
    
    int nx = 4 * Comm.NumProc();
    int ny = 4; 
    int nz = 4;
    int mx = Comm.NumProc(), my = 1, mz = 1;

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
                BoundaryValue, BoundaryType);

    LinearProblem FiniteElementProblem(Grid, Laplacian, A, LHS, RHS); 
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
    
    // ========================= //
    // After the solution:       //
    // - computation of the norm //
    // - output using MEDIT      //
    // ========================= //
    
    FiniteElementProblem.ComputeNorms(LHS, ExactSolution);
    
    MEDITInterface MEDIT(Comm);
    MEDIT.Write(Grid, "Laplacian3D", LHS);

  }
  catch (Exception& rhs)
  {
    if (Comm.MyPID() == 0)
      rhs.Print();
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
