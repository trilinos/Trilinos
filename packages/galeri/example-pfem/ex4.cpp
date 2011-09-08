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

#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsGraph.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FEVector.h"

#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Generator.h"
#include "Galeri_grid_Rebalance.h"
#include "Galeri_quadrature_Segment.h"
#include "Galeri_problem_VectorLaplacian.h"
#include "Galeri_viz_MEDIT.h"

using namespace Galeri;

class MyVectorLaplacian 
{
  public:
    static inline double 
    getElementLHS(const double& x, const double& y, const double& z,
                  const int& ieq, const int& jeq,
                  const double& phi_i,
                  const double& phi_i_derx, 
                  const double& phi_i_dery,
                  const double& phi_i_derz,
                  const double& phi_j,
                  const double& phi_j_derx,
                  const double& phi_j_dery,
                  const double& phi_j_derz)
    {
      if (ieq == 0 &&  jeq == 0)
        return(getEpsilon(x, y, z, 0) * (phi_i_derx * phi_j_derx + 
                                      phi_i_dery * phi_j_dery + 
                                      phi_i_derz * phi_j_derz));
      else if (ieq == 0 && jeq == 1)
        return(phi_i * phi_j);
      else if (ieq == 1 &&  jeq == 1)
        return(getEpsilon(x, y, z, 1) * (phi_i_derx * phi_j_derx + 
                                      phi_i_dery * phi_j_dery + 
                                      phi_i_derz * phi_j_derz));
      else if (ieq == 1 && jeq == 0)
        return(phi_i * phi_j);
      else
        return(0.0);
    }

    static inline double
    getElementRHS(const double& x, const double& y, const double& z,
                  const int &eq, const double& phi_i)

    {
      if (eq == 0) 
        return((-2.0 * getEpsilon(x, y, z, 0) * exp(x) * exp(y) +
                getExactSolution('f', x, y, z, 1)) * phi_i);
      else
        return(getExactSolution('f', x, y, z, 0) * phi_i);
    }

    static inline double
    getBoundaryValue(const double& x, const double& y, const double& z,
                     const int& eq)
    {
      return(getExactSolution('f', x, y, z, eq));
    }

    static inline char
    getBoundaryType(const int ID, const double& x, const double& y, const double& z,
                    const int& eq)
    {
      return('d');
    }

    static inline double 
    getExactSolution(const char& what, const double& x, 
                     const double& y, const double& z, const int eq = 0)
    {
      if (eq == 0)
      {
        if (what == 'f' || what == 'x' || what =='y')
          return(exp(x) * exp(y));
        else
          return (0.0);
      }
      else if (eq == 1)
      {
        if (what == 'f')
          return(x + y);
        else if (what == 'x' || what == 'y')
          return(1.0);
        else
          return (0.0);
      }
      else
        return(0.0);
    }

    static inline double
    getEpsilon(const double& x, const double& y, const double& z, const int& eq)
    {
      if (eq == 0) return(1.0);
      else         return(10.);
    }
};

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  


  Galeri::core::Workspace::setNumDimensions(2);

  Galeri::grid::Loadable domain, boundary;

  int numGlobalElementsX = 64 * comm.NumProc();
  int numGlobalElementsY = 64;

  int mx = comm.NumProc();
  int my = 1;

  Galeri::grid::Generator::
  getSquareWithTriangles(comm, numGlobalElementsX, numGlobalElementsY,
                    mx, my, domain, boundary);

  // ============================================================ //
  // We are now ready to create the linear problem.               //
  // First, we need to define the Epetra_Map for the matrix,      //
  // where each grid vertex is assigned to a different            //
  // processor. To keep things simple, we use a linear partition. //
  // Then, we allocate the matrix (A), the solution vector (LHS), //
  // and the right-hand side (RHS).                               //
  // ============================================================ //
  
  int numPDEs = 1;

  Galeri::problem::VectorLaplacian<MyVectorLaplacian> problem(numPDEs, "Triangle");

  Epetra_BlockMap linearMatrixMap(domain.getNumGlobalVertices(), numPDEs, 0, comm);

  Epetra_FECrsGraph linearGraph(Copy, linearMatrixMap, 0);

  problem.createGraph(domain, linearGraph);

  Epetra_FECrsGraph* graph = Galeri::grid::Rebalance::balanceGraph(linearGraph);

  Epetra_FEVbrMatrix A(Copy, *graph);
  Epetra_FEVector    LHS(graph->Map());
  Epetra_FEVector    RHS(graph->Map());

  problem.integrate(domain, A, RHS);

  LHS.PutScalar(0.0);

  problem.imposeDirichletBoundaryConditions(boundary, A, RHS, LHS);

#if 0
  // ============================================================ //
  // Solving the linear system is the next step, quite easy       //
  // because we just call AztecOO and we wait for the solution... //
  // ============================================================ //
  
  Epetra_LinearProblem linearProblem(&A, &LHS, &RHS);
  AztecOO solver(linearProblem);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);
  solver.SetAztecOption(AZ_output, 16);

  solver.Iterate(1550, 1e-9);
#endif

  Epetra_MultiVector* LHSComponent = Galeri::core::Workspace::createMultiVectorComponent(LHS);

  for (int ieq = 0; ieq < numPDEs; ++ieq)
  {
    Galeri::core::Workspace::extractMultiVectorComponent(LHS, ieq, *LHSComponent);

    char fileName[80];
    sprintf(fileName, "sol%d", ieq);
    Galeri::viz::MEDIT::Write(comm, domain, fileName, *LHSComponent);

    problem.setEquation(ieq);
    problem.computeNorms(domain, *LHSComponent);
  }
  

  delete LHSComponent;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
