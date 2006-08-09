// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Segment.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Rebalance.h"
#include "Galeri_quadrature_Segment.h"
#include "Galeri_problem_ScalarLaplacian.h"

using namespace Galeri;

class MyScalarLaplacian 
{
  public:
    static inline double 
    getElementLHS(const double& x, 
                  const double& y, 
                  const double& z,
                  const double& phi_i,
                  const double& phi_i_derx, 
                  const double& phi_i_dery,
                  const double& phi_i_derz,
                  const double& phi_j,
                  const double& phi_j_derx,
                  const double& phi_j_dery,
                  const double& phi_j_derz)
    {
      return(phi_i_derx * phi_j_derx + 
             phi_i_dery * phi_j_dery + 
             phi_i_derz * phi_j_derz);
    }

    static inline double
    getElementRHS(const double& x,
                  const double& y,
                  const double& z,
                  const double& phi_i)

    {
      return(2.0 * phi_i);
    }

    static inline double
    getBoundaryValue(const double& x, const double& y, const double& z)
    {
      return(0.0);
    }

    static inline char
    getBoundaryType(const int ID, const double& x, const double& y, const double& z)
    {
      return('d');
    }

    static inline double 
    getExactSolution(const char& what, const double& x, 
                     const double& y, const double& z)
    {
      if (what == 'f')
        return(x * (1 - x));
      else if (what == 'x')
        return(1.0 - 2.0 * x);
      else
        return(0.0);
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

  // ======================================================= //
  // Specifies the dimensionality of the problem: 1, 2, or 3 //
  // Creates a 1D grid on (0, 1) composed by segments. Each  //
  // processor will have 4 elements.                         //
  // ======================================================= //
  
  Galeri::core::Workspace::setNumDimensions(1);

  int numMyElements = 4;
  int numGlobalElements = numMyElements * comm.NumProc();

  Galeri::grid::Loadable domain(comm, numGlobalElements, 
                             numMyElements, "Segment");

  // ===================================================== //
  // Each processor inserts locally owned elements, then   //
  // call freezeCoordinates(), which computes the set      //
  // of locally owned vertices. Then, the coordinates of   //
  // these vertices is inserted. Finanlly, the grid is     //
  // freezed by calling freezeConnectivity().              //
  // ===================================================== //
  
  for (int LEID = 0; LEID < numMyElements; ++LEID)
  {
    int GEID = domain.getGEID(LEID);

    domain.setGlobalConnectivity(GEID, 0, GEID);
    domain.setGlobalConnectivity(GEID, 1, GEID + 1);
  }

  domain.freezeConnectivity();

  double h = 1.0 / numGlobalElements;

  for (int LVID = 0; LVID < domain.getNumMyVertices(); ++LVID)
  {
    int GVID = domain.getGVID(LVID);
    domain.setGlobalCoordinates(GVID, 0, h * GVID);
  }

  domain.freezeCoordinates();

  // ========================================================= //
  // We now create the set of boundary nodes. For simplicity,  //
  // both nodes (the one on the left and the one on the right) //
  // are defined on processor 0. The nodes have coordinates    //
  // (0.0) and (1.0), and they are of Dirichlet type.          //
  // ========================================================= //

  int numMyBoundaryElements = (comm.MyPID() == 0)?2:0; 

  Galeri::grid::Loadable boundary(comm, 2, numMyBoundaryElements, "Point");

  if (comm.MyPID() == 0)
  {
    boundary.setGlobalConnectivity(0, 0, 0);
    boundary.setGlobalConnectivity(1, 0, domain.getNumGlobalElements());
  }

  boundary.freezeConnectivity();

  if (comm.MyPID() == 0)
  {
    boundary.setGlobalCoordinates(0, 0, 0.0);
    boundary.setGlobalCoordinates(1, 0, 1.0);
  }

  boundary.freezeCoordinates();

  // ============================================================ //
  // We are now ready to create the linear problem.               //
  // First, we need to define the Epetra_Map for the matrix,      //
  // where each grid vertex is assigned to a different            //
  // processor. To keep things simple, we use a linear partition. //
  // Then, we allocate the matrix (A), the solution vector (LHS), //
  // and the right-hand side (RHS).                               //
  // ============================================================ //
  
  Epetra_Map matrixMap(domain.getNumGlobalElements() + 1, 0, comm);

  Epetra_FECrsMatrix A(Copy, matrixMap, 0);
  Epetra_FEVector    LHS(matrixMap);
  Epetra_FEVector    RHS(matrixMap);

  Galeri::problem::ScalarLaplacian<MyScalarLaplacian> problem("Segment", 1, 4);

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

  solver.Iterate(150, 1e-9);
#endif

  // now compute the norm of the solution
  
  problem.computeNorms(domain, LHS);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
