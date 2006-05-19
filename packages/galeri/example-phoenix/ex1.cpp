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
#include "../src-phoenix/phx_grid_Triangle.h"
#include "../src-phoenix/phx_grid_Quad.h"
#include "../src-phoenix/phx_grid_Loadable.h"
#include "../src-phoenix/phx_grid_SerialXML.h"
#include "../src-phoenix/phx_grid_Generator.h"
#include "../src-phoenix/phx_quadrature_Element.h"
#include "../src-phoenix/phx_quadrature_Triangle.h"
#include "../src-phoenix/phx_quadrature_Quad.h"
#include "../src-phoenix/phx_problem_ScalarLaplacian.h"
#include "../src-phoenix/phx_viz_MEDIT.h"
#include "../src-phoenix/phx_viz_VTK.h"
#include <fstream>

using namespace phx;

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

  string XMLFileName = "/Users/marzio/test.xml";

  map<string, RefCountPtr<phx::grid::Loadable> > patches;
  //patches = phx::grid::SerialXML::readFile(comm, XMLFileName);
  patches = phx::grid::Generator::getSquareWithQuad(comm, 10, 10, 1, 1);

  RefCountPtr<phx::grid::Loadable> domain          = patches["domain"];
  RefCountPtr<phx::grid::Loadable> boundary        = patches["boundary"];
  RefCountPtr<phx::grid::Element>  domainElement   = domain->getElement();
  RefCountPtr<phx::grid::Element>  boundaryElement = boundary->getElement();

  RefCountPtr<Epetra_Map> matrixMap = phx::core::Utils::createMatrixMap(comm, patches);

  Epetra_FECrsMatrix matrix(Copy, *matrixMap, 0);
  Epetra_FEVector    lhs(*matrixMap);
  Epetra_FEVector    rhs(*matrixMap);

  Epetra_IntSerialDenseVector vertexList(domainElement->getNumVertices());

  phx::quadrature::Quad domainQuadrature(1);
  phx::problem::ScalarLaplacian problem;
  Epetra_SerialDenseMatrix elementLHS(domainElement->getNumVertices(),domainElement->getNumVertices());
  Epetra_SerialDenseVector elementRHS(domainElement->getNumVertices());

  for (int i = 0; i < domain->getNumMyElements(); ++i)
  {
    // load the element vertex IDs
    for (int j = 0; j < domainElement->getNumVertices(); ++j)
      vertexList[j] = domain->getMyConnectivity(i, j);

    // load the element coordinates
    for (int j = 0; j < domainElement->getNumVertices(); ++j)
      for (int k = 0; k < domainElement->getNumDimensions(); ++k) 
        domainQuadrature(j, k) = domain->getGlobalCoordinates(vertexList[j], k);

    problem.Integrate(domainQuadrature, elementLHS, elementRHS);

    matrix.InsertGlobalValues(vertexList, elementLHS);
    rhs.SumIntoGlobalValues(vertexList, elementRHS);
  }

  matrix.GlobalAssemble();
  rhs.GlobalAssemble();

  Teuchos::Hashtable<int, double> dirichletRows;

  vector<double> coord(boundaryElement->getNumDimensions());
  
  for (int i = 0; i < boundary->getNumMyElements(); ++i)
  {
    for (int j = 0; j < boundaryElement->getNumVertices(); ++j)
    {
      int GID = boundary->getGlobalConnectivity(i, j);
      for (int k = 0; k < boundaryElement->getNumDimensions(); ++k)
        coord[k] = boundary->getGlobalCoordinates(GID, k);
      
      dirichletRows.put(GID, 0.0);
    }
  }

  for (int i = 0; i < matrix.NumMyRows(); ++i)
  {
    int GID = matrix.RowMatrixRowMap().GID(i);

    bool isDirichlet = false;
    if (dirichletRows.containsKey(GID)) isDirichlet = true;

    int* indices;
    double* values;
    int numEntries;
    matrix.ExtractMyRowView(i, numEntries, values, indices);

    if (isDirichlet)
    {
      for (int j = 0; j < numEntries; ++j)
      if (indices[j] != i) values[j] = 0.0;
      else values[j] = 1.0;
      rhs[0][i] = dirichletRows.get(GID);
    }
    else
    {
      for (int j = 0; j < numEntries; ++j)
      {
        if (indices[j] == i) continue;
        if (dirichletRows.containsKey(matrix.RowMatrixColMap().GID(indices[j]))) values[j] = 0.0;
        // FIXME: RHS??
      }
    }
  }

  lhs.PutScalar(0.0);

  AztecOO solver(&matrix, &lhs, &rhs);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);
  solver.SetAztecOption(AZ_output, 16);

  solver.Iterate(150, 1e-9);

  phx::viz::VTK::Write(comm, patches["domain"], "sol", lhs);

  Epetra_SerialDenseVector elementSol(domainElement->getNumVertices());
  Epetra_SerialDenseVector elementNorm(2);
  double normL2 = 0.0, semiNormH1 = 0.0;

  for (int i = 0; i < domain->getNumMyElements(); ++i)
  {
    for (int j = 0; j < domainElement->getNumVertices(); ++j)
    {
      vertexList[j] = domain->getMyConnectivity(i, j);
      elementSol[j] = lhs[0][vertexList[j]];
    }

    // load the element coordinates
    for (int j = 0; j < domainElement->getNumVertices(); ++j)
      for (int k = 0; k < domainElement->getNumDimensions(); ++k) 
        domainQuadrature(j, k) = domain->getGlobalCoordinates(vertexList[j], k);

    problem.computeNorm(domainQuadrature, elementSol, elementNorm);

    normL2     += elementNorm[0];
    semiNormH1 += elementNorm[1];
  }


  cout << "Norm L2 = " << normL2 << endl;
  cout << "SemiNorm H1 = " << semiNormH1 << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
