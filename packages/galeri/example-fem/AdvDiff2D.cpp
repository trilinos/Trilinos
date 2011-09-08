// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
//   - \mu \nabla u + c_x * u_x + c_y * u_y = f    on \Omega
//                                        u = g    on \partial \Omega
// 
// where \Omega is a 2D rectangle, divided into triangles.
// `f' is specified by function `Force()', the Dirichlet boundary condition
// by function `BoundaryValue()', and the value of \mu and
// c_x and c_y can be changed in the functions Diffusion() and
// ConvX() and ConvY(). The code solves the corresponding 
// linear system using a simple LAPACK interface, and writes
// the solution inot a MEDIT-compatible format.
//
// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 15-Sep-05.
// ==========================================================

double Diffusion(const double& x, const double& y, const double& z)
{
  return (1.0);
}

double conv = 5000;
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

int BoundaryType(const int& Patch)
{
  return(GALERI_DIRICHLET);
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
    FileGrid Grid(Comm, "Square.grid");

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

    SUPGVariational<TriangleQuadrature>
      AdvDiff(NumQuadratureNodes, Diffusion, ConvX, ConvY, ConvZ, 
              Source, Force, BoundaryValue, BoundaryType);

    LinearProblem FiniteElementProblem(Grid, AdvDiff, A, LHS, RHS); 
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
    MEDIT.Write(Grid, "AdvDiff2D", LHS);

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
