// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <limits>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Galeri_XpetraMaps.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

using Teuchos::RCP;

/////////////////////////////////////////////////////////////////////////////
// Program to demonstrate use of Zoltan2 to partition a TPetra matrix
// using graph partitioning via Scotch or ParMETIS.
/////////////////////////////////////////////////////////////////////////////

int main(int narg, char** arg)
{
  // Establish session; works both for MPI and non-MPI builds
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  // Useful typedefs:  Tpetra types
  // In this example, we'll use Tpetra defaults for local/global ID type
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef Tpetra::Details::DefaultTypes::scalar_type scalar_t;
  typedef Tpetra::CrsMatrix<scalar_t, localId_t, globalId_t> Matrix_t;
  typedef Tpetra::MultiVector<scalar_t, localId_t, globalId_t> MultiVector_t;
  typedef Tpetra::Vector<scalar_t, localId_t, globalId_t> Vector_t;

  // Useful typedefs:  Zoltan2 types
  typedef Zoltan2::XpetraCrsMatrixAdapter<Matrix_t> MatrixAdapter_t;
  typedef Zoltan2::XpetraMultiVectorAdapter<Vector_t> MultiVectorAdapter_t;

  // Input parameters with default values
  std::string method = "scotch";    // Partitioning method
  globalId_t nx = 50, ny = 40, nz = 30; // Dimensions of mesh corresponding to
                                    // the matrix to be partitioned

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("method", &method,
                 "Partitioning method to use:  scotch or parmetis.");
  cmdp.setOption("nx", &nx,
                 "number of gridpoints in X dimension for "
                 "mesh used to generate matrix; must be >= 1.");
  cmdp.setOption("ny", &ny,
                 "number of gridpoints in Y dimension for "
                 "mesh used to generate matrix; must be >= 1.");
  cmdp.setOption("nz", &nz,
                 "number of gridpoints in Z dimension for "
                 "mesh used to generate matrix; must be >= 1.");
  cmdp.parse(narg, arg);

  if ((nx < 1) || (ny < 1) || (nz < 1)) {
    std::cout << "Input error:  nx, ny and nz must be >= 1" << std::endl;
    return -1;
  }

  // For this example, generate a matrix using Galeri with the
  // default Tpetra distribution:
  //   Laplace3D matrix corresponding to mesh with dimensions (nx X ny X nz)
  //   with block row-based distribution
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("nz", nz);
  Tpetra::global_size_t nGlobalElements = nx * ny * nz;

  RCP<Matrix_t> origMatrix;

  try {
    RCP<const Map_t> map = rcp(new Map_t(nGlobalElements, 0, comm));

    typedef Galeri::Xpetra::Problem<Map_t,Matrix_t,MultiVector_t> Galeri_t;
    RCP<Galeri_t> galeriProblem =
                  Galeri::Xpetra::BuildProblem<scalar_t, localId_t, globalId_t,
                                     Map_t, Matrix_t, MultiVector_t>
                                     ("Laplace3D", map, galeriList);
    origMatrix = galeriProblem->BuildMatrix();
  }
  catch (std::exception &e) {
    std::cout << "Exception in Galeri matrix generation. " << e.what() << std::endl;
    return -1;
  }

  if (me == 0)
    std::cout << "NumRows     = " << origMatrix->getGlobalNumRows() << std::endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << std::endl
         << "NumProcs = " << comm->getSize() << std::endl;

  // Create vectors to use with the matrix for sparse matvec.
  RCP<Vector_t> origVector, origProd;
  origProd   = Tpetra::createVector<scalar_t,localId_t,globalId_t>(
                                    origMatrix->getRangeMap());
  origVector = Tpetra::createVector<scalar_t,localId_t,globalId_t>(
                                    origMatrix->getDomainMap());
  origVector->randomize();

  // Specify partitioning parameters
  Teuchos::ParameterList param;
  param.set("partitioning_approach", "partition");
  param.set("algorithm", method);

  // Create an input adapter for the Tpetra matrix.
  MatrixAdapter_t adapter(origMatrix);

  // Create and solve partitioning problem
  Zoltan2::PartitioningProblem<MatrixAdapter_t> problem(&adapter, &param);

  try {
    problem.solve();
  }
  catch (std::exception &e) {
    std::cout << "Exception returned from solve(). " << e.what() << std::endl;
    return -1;
  }

  // Redistribute matrix and vector into new matrix and vector.
  // Can use PartitioningSolution from matrix to redistribute the vectors, too.

  if (me == 0) std::cout << "Redistributing matrix..." << std::endl;
  RCP<Matrix_t> redistribMatrix;
  adapter.applyPartitioningSolution(*origMatrix, redistribMatrix,
                                    problem.getSolution());

  if (me == 0) std::cout << "Redistributing vectors..." << std::endl;
  RCP<Vector_t> redistribVector;
  MultiVectorAdapter_t adapterVector(origVector);
  adapterVector.applyPartitioningSolution(*origVector, redistribVector,
                                          problem.getSolution());

  // Create a new product vector for sparse matvec
  RCP<Vector_t> redistribProd;
  redistribProd = Tpetra::createVector<scalar_t,localId_t,globalId_t>(
                                       redistribMatrix->getRangeMap());

  // SANITY CHECK
  // A little output for small problems
  if (origMatrix->getGlobalNumRows() <= 50) {
    std::cout << me << " ORIGINAL:  ";
    for (size_t i = 0; i < origVector->getLocalLength(); i++)
      std::cout << origVector->getMap()->getGlobalElement(i) << " ";
    std::cout << std::endl;
    std::cout << me << " REDISTRIB: ";
    for (size_t i = 0; i < redistribVector->getLocalLength(); i++)
      std::cout << redistribVector->getMap()->getGlobalElement(i) << " ";
    std::cout << std::endl;
  }

  // SANITY CHECK
  // check that redistribution is "correct"; perform matvec with
  // original and redistributed matrices/vectors and compare norms.

  if (me == 0) std::cout << "Matvec original..." << std::endl;
  origMatrix->apply(*origVector, *origProd);
  scalar_t origNorm = origProd->norm2();
  if (me == 0)
    std::cout << "Norm of Original matvec prod:       " << origNorm << std::endl;

  if (me == 0) std::cout << "Matvec redistributed..." << std::endl;
  redistribMatrix->apply(*redistribVector, *redistribProd);
  scalar_t redistribNorm = redistribProd->norm2();
  if (me == 0)
    std::cout << "Norm of Redistributed matvec prod:  " << redistribNorm << std::endl;

  if (me == 0) {
    const scalar_t epsilon = 0.001;
    if (redistribNorm > origNorm+epsilon || redistribNorm < origNorm-epsilon)
      std::cout << "Mat-Vec product changed; FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return 0;
}
