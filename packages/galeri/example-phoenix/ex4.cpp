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
#include "AztecOO.h"

#include "phx_core_Object.h"
#include "phx_core_Utils.h"
#include "phx_grid_Triangle.h"
#include "phx_grid_Loadable.h"
#include "phx_grid_Generator.h"
#include "phx_quadrature_Segment.h"
#include "phx_problem_VectorLaplacian.h"
#include "phx_viz_MEDIT.h"

using namespace phx;

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
  


  phx::core::Utils::setNumDimensions(2);

  phx::grid::Loadable domain, boundary;

  int numGlobalElementsX = 64 * comm.NumProc();
  int numGlobalElementsY = 64;

  int mx = comm.NumProc();
  int my = 1;

  phx::grid::Generator::
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
  
  int numPDEs = 2;

  phx::problem::VectorLaplacian<MyVectorLaplacian> problem(numPDEs, "Triangle");

  Epetra_BlockMap matrixMap(domain.getNumGlobalElements() + 1, numPDEs, 0, comm);

  Epetra_FECrsGraph  Graph(Copy, matrixMap, 0);

  problem.createGraph(domain, Graph);

  Epetra_FEVbrMatrix A(Copy, Graph);
  Epetra_FEVector    LHS(matrixMap);
  Epetra_FEVector    RHS(matrixMap);

  problem.integrate(domain, A, RHS);

  LHS.PutScalar(0.0);

  problem.imposeDirichletBoundaryConditions(boundary, A, RHS, LHS);

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

  Epetra_MultiVector* LHSComponent = phx::core::Utils::createMultiVectorComponent(LHS);

  for (int ieq = 0; ieq < numPDEs; ++ieq)
  {
    phx::core::Utils::extractMultiVectorComponent(LHS, ieq, *LHSComponent);

    char fileName[80];
    sprintf(fileName, "sol%d", ieq);
    phx::viz::MEDIT::Write(comm, domain, fileName, *LHSComponent);

    problem.setEquation(ieq);
    problem.computeNorms(domain, *LHSComponent);
  }
  

  delete LHSComponent;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
