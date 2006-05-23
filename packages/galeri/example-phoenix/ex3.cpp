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

  phx::core::Utils::setNumDimensions(3);

  phx::grid::Loadable domain, boundary;

  map<char, int> numGlobalElements;
  numGlobalElements['x'] = 20;
  numGlobalElements['y'] = 20;
  numGlobalElements['z'] = 20;
  numGlobalElements['p'] = numGlobalElements['x'] * numGlobalElements['y'];
  numGlobalElements['q'] = numGlobalElements['y'] * numGlobalElements['z'];
  numGlobalElements['r'] = numGlobalElements['x'] * numGlobalElements['z'];
  numGlobalElements['a'] = numGlobalElements['x'] * numGlobalElements['y'] * numGlobalElements['z'];

  map<char, int> numGlobalVertices;
  numGlobalVertices['x'] = numGlobalElements['x'] + 1;
  numGlobalVertices['y'] = numGlobalElements['y'] + 1;
  numGlobalVertices['z'] = numGlobalElements['z'] + 1;
  numGlobalVertices['p'] = numGlobalVertices['x'] * numGlobalVertices['y'];
  numGlobalVertices['q'] = numGlobalVertices['y'] * numGlobalVertices['z'];
  numGlobalVertices['r'] = numGlobalVertices['x'] * numGlobalVertices['z'];
  numGlobalVertices['a'] = numGlobalVertices['x'] * numGlobalVertices['y'] * numGlobalVertices['z'];

  map<char, double> length;
  length['x'] = 1.0;
  length['y'] = 1.0;
  length['z'] = 1.0;

  domain.initialize(comm, numGlobalElements['a'], -1, "Hex");

  for (int iz = 0; iz < numGlobalElements['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalElements['y']; ++iy)
    {
      for (int ix = 0; ix < numGlobalElements['x']; ++ix)
      {
        int GEID = iz * numGlobalElements['p'] + iy * numGlobalElements['x'] + ix;
        int offset = iz * numGlobalVertices['p'] + iy * numGlobalVertices['x'] + ix;

        domain.setGlobalConnectivity(GEID, 0, offset);
        domain.setGlobalConnectivity(GEID, 1, offset + 1);
        domain.setGlobalConnectivity(GEID, 2, offset + numGlobalVertices['x'] + 1);
        domain.setGlobalConnectivity(GEID, 3, offset + numGlobalVertices['x']);
        domain.setGlobalConnectivity(GEID, 4, offset + numGlobalVertices['p']);
        domain.setGlobalConnectivity(GEID, 5, offset + numGlobalVertices['p'] + 1);
        domain.setGlobalConnectivity(GEID, 6, offset + numGlobalVertices['p'] + numGlobalVertices['x'] + 1);
        domain.setGlobalConnectivity(GEID, 7, offset + numGlobalVertices['p'] + numGlobalVertices['x']);
      }
    }
  }

  domain.freezeConnectivity();

  double hx = length['x'] / numGlobalElements['x'];
  double hy = length['y'] / numGlobalElements['y'];
  double hz = length['z'] / numGlobalElements['z'];

  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
      {
        int GVID = iz * numGlobalVertices['p'] + iy * numGlobalVertices['x'] + ix;
        domain.setGlobalCoordinates(GVID, 0, hx * ix);
        domain.setGlobalCoordinates(GVID, 1, hy * iy);
        domain.setGlobalCoordinates(GVID, 2, hz * iz);
      }
    }
  }

  domain.freezeCoordinates();

  int numGlobalBoundaries = 2 * numGlobalVertices['p'] +
                            2 * numGlobalVertices['q'] +
                            2 * numGlobalVertices['r'];
  boundary.initialize(comm, numGlobalBoundaries, -1, "Point");

  int count = 0;

  // bottom 
  for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iy * numGlobalVertices['x'] + ix;
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  // top 
  for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iy * numGlobalVertices['x'] + ix +
        numGlobalVertices['p'] * numGlobalElements['z'];
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  // front 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iz * numGlobalVertices['p'] + ix;
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  // rear 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = (iz + 1) * numGlobalVertices['p'] - numGlobalVertices['x'] + ix;
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  // left 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
    {
      int GBID = iz * numGlobalVertices['p'] + iy * numGlobalVertices['x'];
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  // right 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
    {
      int GBID = iz * numGlobalVertices['p'] + (iy + 1) * numGlobalVertices['x'] - 1;
      boundary.setGlobalConnectivity(count++, 0, GBID);
    }
  }

  boundary.freezeConnectivity();

  // bottom
  for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iy * numGlobalVertices['x'] + ix;

      boundary.setGlobalCoordinates(GBID, 0, hx * ix);
      boundary.setGlobalCoordinates(GBID, 1, hy * iy);
      boundary.setGlobalCoordinates(GBID, 2, 0.0);
    }
  }

  // top
  for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iy * numGlobalVertices['x'] + ix +
        numGlobalVertices['p'] * numGlobalElements['z'];

      boundary.setGlobalCoordinates(GBID, 0, hx * ix);
      boundary.setGlobalCoordinates(GBID, 1, hy * iy);
      boundary.setGlobalCoordinates(GBID, 2, length['z']);
    }
  }

  // front
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = iz * numGlobalVertices['p'] + ix;

      boundary.setGlobalCoordinates(GBID, 0, hx * ix);
      boundary.setGlobalCoordinates(GBID, 1, 0.0);
      boundary.setGlobalCoordinates(GBID, 2, hz * iz);
    }
  }

  // rear
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numGlobalVertices['x']; ++ix)
    {
      int GBID = (iz + 1) * numGlobalVertices['p'] - numGlobalVertices['x'] + ix;

      boundary.setGlobalCoordinates(GBID, 0, hx * ix);
      boundary.setGlobalCoordinates(GBID, 1, length['y']);
      boundary.setGlobalCoordinates(GBID, 2, hz * iz);
    }
  }

  // left 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
    {
      int GBID = iz * numGlobalVertices['p'] + iy * numGlobalVertices['x'];

      boundary.setGlobalCoordinates(GBID, 0, 0.0);
      boundary.setGlobalCoordinates(GBID, 1, hy * iy);
      boundary.setGlobalCoordinates(GBID, 2, hz * iz);
    }
  }

  // right 
  for (int iz = 0; iz < numGlobalVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numGlobalVertices['y']; ++iy)
    {
      int GBID = iz * numGlobalVertices['p'] + (iy + 1) * numGlobalVertices['x'] - 1;

      boundary.setGlobalCoordinates(GBID, 0, length['x']);
      boundary.setGlobalCoordinates(GBID, 1, hy * iy);
      boundary.setGlobalCoordinates(GBID, 2, hz * iz);
    }
  }

  boundary.freezeCoordinates();

  Epetra_Map matrixMap(domain.getNumGlobalVertices(), 0, comm);

  Epetra_FECrsMatrix A(Copy, matrixMap, 0);
  Epetra_FEVector    LHS(matrixMap);
  Epetra_FEVector    RHS(matrixMap);

  phx::problem::ScalarLaplacian<Laplacian> problem("Hex", 1, 8);

  problem.integrate(domain, A, RHS);

  LHS.PutScalar(0.0);

  problem.imposeDirichletBoundaryConditions(boundary, A, RHS, LHS);

  // ============================================================ //
  // Solving the linear system is the next step, quite easy       //
  // because we just call AztecOO and we wait for the solution... //
  // ============================================================ //
  
  Epetra_LinearProblem linearProblem(&A, &LHS, &RHS);
  AztecOO solver(linearProblem);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  solver.SetAztecOption(AZ_output, 16);

  solver.Iterate(1550, 1e-9);

  phx::viz::MEDIT::Write(comm, domain, "sol", LHS);

  // now compute the norm of the solution
  
  problem.computeNorms(domain, LHS);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
