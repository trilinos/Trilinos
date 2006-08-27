// @HEADER
// ***********************************************************************
// 
//                IFPACK: Robust Algebraic Preconditioning Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
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
#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Generator.h"
#include "Galeri_quadrature_Hex.h"
#include "Galeri_problem_ScalarLaplacian.h"
#include "Galeri_viz_MEDIT.h"

#include "AztecOO.h"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

using namespace Galeri;

class Laplacian 
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
      return(-getExactSolution('f', x, y, z) * phi_i);
    }

    static inline double
    getBoundaryValue(const double& x, const double& y, const double& z)
    {
      return(getExactSolution('f', x, y, z));
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
        return(exp(x) + exp(y) + exp(z));
      else if (what == 'x')
        return(exp(x));
      else if (what == 'y')
        return(exp(y));
      else if (what == 'z')
        return(exp(z));
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

  Galeri::core::Workspace::setNumDimensions(3);

  Galeri::grid::Loadable domain, boundary;

  int numGlobalElementsX = 2 * comm.NumProc();
  int numGlobalElementsY = 2;
  int numGlobalElementsZ = 2;

  int mx = comm.NumProc();
  int my = 1;
  int mz = 1;

  Galeri::grid::Generator::
  getCubeWithHexs(comm, numGlobalElementsX, numGlobalElementsY, numGlobalElementsZ,
                  mx, my, mz, domain, boundary);

  Epetra_Map matrixMap(domain.getNumGlobalVertices(), 0, comm);

  Epetra_FECrsMatrix A(Copy, matrixMap, 0);
  Epetra_FEVector    LHS(matrixMap);
  Epetra_FEVector    RHS(matrixMap);

  Galeri::problem::ScalarLaplacian<Laplacian> problem("Hex", 1, 8);

  problem.integrate(domain, A, RHS);

  LHS.PutScalar(0.0);

  problem.imposeDirichletBoundaryConditions(boundary, A, RHS, LHS);

  // ============================================================ //
  // Solving the linear system is the next step, using the IFPACK //
  // factory. This is done by using the IFPACK factory, then      //
  // asking for IC preconditioner, and setting few parameters     //
  // using a Teuchos::ParameterList.                              //
  // ============================================================ //
  
  Ifpack Factory;
  Ifpack_Preconditioner* Prec = Factory.Create("IC", &A, 0);

  Teuchos::ParameterList list;
  
  list.set("fact: level-of-fill", 1);
  IFPACK_CHK_ERR(Prec->SetParameters(list));
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());

  Epetra_LinearProblem linearProblem(&A, &LHS, &RHS);

  AztecOO solver(linearProblem);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetPrecOperator(Prec);
  solver.Iterate(1550, 1e-9);

  // visualization using MEDIT -- a VTK module is avaiable as well
  Galeri::viz::MEDIT::write(domain, "sol", LHS);

  // now compute the norm of the solution
  problem.computeNorms(domain, LHS);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
