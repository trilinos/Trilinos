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
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"

#include "Galeri_core_Constants.h"
#include "Galeri_core_Object.h"
#include "Galeri_core_Utils.h"
#include "Galeri_grid_Segment.h"
#include "Galeri_grid_Quad.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Generator.h"
#include "Galeri_quadrature_Quad.h"
#include "Galeri_problem_ScalarLaplacian.h"
#include "Galeri_viz_MEDIT.h"

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
      return(2.0 * (y * (1 - y) + x * (1 - x)) * phi_i);
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
        return(x * (1 - x) * y * (1 - y));
      else if (what == 'x')
        return((1.0 -2.0 * x) * y * (1 - y));
      else if (what == 'y')
        return(x * (1 - x) * (1.0 - 2.0 * y));
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

  Galeri::core::Utils::setNumDimensions(2);

  Galeri::grid::Loadable domain, boundary;

  int numGlobalElementsX = 4 * comm.NumProc();
  int numGlobalElementsY = numGlobalElementsX;

  int mx = comm.NumProc();
  int my = 1;

#if 1
  Galeri::grid::Generator::
  getSquareWithTriangles(comm, numGlobalElementsX, numGlobalElementsY,
                    mx, my, domain, boundary);
#else
  string XMLFileName = "/Users/marzio/test.xml";


  Galeri::grid::SerialXML XMLReader;
  map<string, Galeri::grid::Loadable> patches = XMLReader.read(comm, XMLFileName);

  domain = patches["domain"];
  boundary = patches["boundary"];
#endif

  Epetra_Map matrixMap(domain.getNumGlobalVertices(), 0, comm);

  Epetra_FECrsMatrix A(Copy, matrixMap, 0);
  Epetra_FEVector    LHS(matrixMap);
  Epetra_FEVector    RHS(matrixMap);

  Galeri::problem::ScalarLaplacian<Laplacian> problem("Triangle");

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
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  solver.SetAztecOption(AZ_output, 16);

  solver.Iterate(1550, 1e-9);
#endif

  Galeri::viz::MEDIT::Write(comm, domain, "sol", LHS);

  // now compute the norm of the solution
  
  problem.computeNorms(domain, LHS);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
