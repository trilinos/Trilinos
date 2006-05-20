#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "AztecOO.h"

#include "../src-phoenix/phx_core_Object.h"
#include "../src-phoenix/phx_core_Utils.h"
#include "../src-phoenix/phx_grid_Element.h"
#include "../src-phoenix/phx_grid_Segment.h"
#include "../src-phoenix/phx_grid_Triangle.h"
#include "../src-phoenix/phx_grid_Quad.h"
#include "../src-phoenix/phx_grid_Tet.h"
#include "../src-phoenix/phx_grid_Hex.h"
#include "../src-phoenix/phx_grid_Loadable.h"
#include "../src-phoenix/phx_grid_SerialXML.h"
#include "../src-phoenix/phx_grid_Generator.h"
#include "../src-phoenix/phx_quadrature_Element.h"
#include "../src-phoenix/phx_quadrature_Segment.h"
#include "../src-phoenix/phx_quadrature_Triangle.h"
#include "../src-phoenix/phx_quadrature_Quad.h"
#include "../src-phoenix/phx_quadrature_Tet.h"
#include "../src-phoenix/phx_quadrature_Hex.h"
#include "../src-phoenix/phx_problem_ScalarLaplacian.h"
#include "../src-phoenix/phx_viz_MEDIT.h"
#include "../src-phoenix/phx_viz_VTK.h"
#include <fstream>

using namespace phx;

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
  // Creates a 1D grid on (0, 1) composed by segments. Each  //
  // processor will have 4 elements. The boundary conditions //
  // are of Dirichlet type.                                  //
  // ======================================================= //
  
  int numDimensions = 1;
  int numMyElements = 32;
  int numGlobalElements = numMyElements * comm.NumProc();

  RefCountPtr<phx::grid::Loadable> domain;
  
  domain = rcp(new phx::grid::Loadable(comm, numGlobalElements, numMyElements, 
                                       "Segment", numDimensions));

  // ===================================================== //
  // Each processor inserts locally owned elements, then   //
  // call freezeCoordinates(), which computes the set      //
  // of locally owned vertices. Then, the coordinates of   //
  // these vertices is inserted. Finanlly, the grid is     //
  // freezed by calling freezeConnectivity().              //
  // ===================================================== //
  
  for (int LEID = 0; LEID < numMyElements; ++LEID)
  {
    int GEID = domain->getGEID(LEID);

    domain->setGlobalConnectivity(GEID, 0, GEID);
    domain->setGlobalConnectivity(GEID, 1, GEID + 1);
  }

  domain->freezeConnectivity();

  double h = 1.0 / numGlobalElements;

  for (int LVID = 0; LVID < domain->getNumMyVertices(); ++LVID)
  {
    int GVID = domain->getGVID(LVID);
    domain->setGlobalCoordinates(GVID, 0, h * GVID);
  }

  domain->freezeCoordinates();

  // The boundary points are now inserted in an equivalent
  // manner.

  int numMyBoundaries = 0;
  if (comm.MyPID() == 0) ++numMyBoundaries; 
  if (comm.MyPID() == comm.NumProc() - 1) ++numMyBoundaries;

  RefCountPtr<phx::grid::Loadable> boundaries;
  boundaries = rcp(new phx::grid::Loadable(comm, 2, numMyBoundaries, "Point", 1));
  if (comm.MyPID() == 0)
    boundaries->setGlobalConnectivity(0, 0, 0);
  if (comm.MyPID() == comm.NumProc() - 1)
    boundaries->setGlobalConnectivity(1, 0, domain->getNumGlobalVertices() - 1);

  boundaries->freezeConnectivity();

  if (comm.MyPID() == 0)
    boundaries->setGlobalCoordinates(0, 0, 0.0);
  if (comm.MyPID() == comm.NumProc() - 1)
    boundaries->setGlobalCoordinates(1, 0, 1.0);

  boundaries->freezeCoordinates();

  // ============================================================ //
  // We are now ready to create the linear problem.               //
  // First, we need to define the Epetra_Map for the matrix,      //
  // where each grid vertex is assigned to a different            //
  // processor. To keep things simple, we use a linear partition. //
  // Then, we allocate the matrix (A), the solution vector (LHS), //
  // and the right-hand side (RHS).                               //
  // ============================================================ //
  
  RefCountPtr<Epetra_Map> matrixMap = 
    rcp(new Epetra_Map(domain->getNumGlobalVertices(), 0, comm));

  phx::quadrature::Segment domainQuadrature(4);
  phx::problem::ScalarLaplacian<MyScalarLaplacian> problem(matrixMap, "Segment", 1, 4);

  problem.integrate(domain, numDimensions);

  problem.imposeDirichletBoundaryConditions(boundaries, numDimensions);

  // ============================================================ //
  // Solving the linear system is the next step, quite easy       //
  // because we just call AztecOO and we wait for the solution... //
  // ============================================================ //
  
  AztecOO solver(problem.getMatrix(), problem.getLHS(), problem.getRHS());
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);
  solver.SetAztecOption(AZ_output, 16);

  solver.Iterate(150, 1e-9);

  // now compute the norm of the solution
  
  problem.computeNorms(domain, numDimensions, *problem.getLHS(), true);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
