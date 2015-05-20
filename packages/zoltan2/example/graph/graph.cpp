// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <limits>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Galeri_XpetraMaps.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

using Teuchos::RCP;
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// Program to demonstrate use of Zoltan2 to partition a TPetra matrix
// using graph partitioning via Scotch or ParMETIS.
/////////////////////////////////////////////////////////////////////////////

int main(int narg, char** arg)
{
  // Establish session; works both for MPI and non-MPI builds
  Teuchos::GlobalMPISession mpiSession(&narg, &arg, NULL);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  int me = comm->getRank();

  // Useful typedefs:  basic types
  typedef int zlno_t;
  typedef long zgno_t;
  typedef double zscalar_t;

  // Useful typedefs:  Tpetra types
  typedef Tpetra::Map<zlno_t, zgno_t> Map_t;
  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> Matrix_t;
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t> MultiVector_t;
  typedef Tpetra::Vector<zscalar_t, zlno_t, zgno_t> Vector_t;

  // Useful typedefs:  Zoltan2 types
  typedef Zoltan2::XpetraCrsMatrixAdapter<Matrix_t> MatrixAdapter_t;
  typedef Zoltan2::XpetraMultiVectorAdapter<Vector_t> MultiVectorAdapter_t;

  // Input parameters with default values
  std::string method = "scotch";    // Partitioning method
  zgno_t nx = 50, ny = 40, nz = 30; // Dimensions of mesh corresponding to
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
    cout << "Input error:  nx, ny and nz must be >= 1" << endl;
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
                  Galeri::Xpetra::BuildProblem<zscalar_t, zlno_t, zgno_t, 
                                     Map_t, Matrix_t, MultiVector_t>
                                     ("Laplace3D", map, galeriList);
    origMatrix = galeriProblem->BuildMatrix();
  }
  catch (std::exception &e) {
    cout << "Exception in Galeri matrix generation. " << e.what() << endl;
    return -1;
  }
  
  if (me == 0)
    cout << "NumRows     = " << origMatrix->getGlobalNumRows() << endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << endl
         << "NumProcs = " << comm->getSize() << endl;

  // Create vectors to use with the matrix for sparse matvec.
  RCP<Vector_t> origVector, origProd;
  origProd   = Tpetra::createVector<zscalar_t,zlno_t,zgno_t>(
                                    origMatrix->getRangeMap());
  origVector = Tpetra::createVector<zscalar_t,zlno_t,zgno_t>(
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
    cout << "Exception returned from solve(). " << e.what() << endl;
    return -1;
  }

  // Basic metric checking of the partitioning solution
  // Not ordinarily done in application code; just doing it for testing here.
  size_t checkNparts = comm->getSize();

  size_t  checkLength = origMatrix->getNodeNumRows();
  const MatrixAdapter_t::part_t *checkParts = 
                                       problem.getSolution().getPartListView();

  // Check for load balance
  size_t *countPerPart = new size_t[checkNparts];
  size_t *globalCountPerPart = new size_t[checkNparts];

  for (size_t i = 0; i < checkNparts; i++) countPerPart[i] = 0;

  for (size_t i = 0; i < checkLength; i++) {
    if (size_t(checkParts[i]) >= checkNparts)
      cout << "Invalid Part " << checkParts[i] << ": FAIL" << endl;
    countPerPart[checkParts[i]]++;
  }
  Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, checkNparts,
                                  countPerPart, globalCountPerPart);

  size_t min = std::numeric_limits<std::size_t>::max();
  size_t max = 0;
  size_t sum = 0;
  size_t minrank = 0, maxrank = 0;
  for (size_t i = 0; i < checkNparts; i++) {
    if (globalCountPerPart[i] < min) {min = globalCountPerPart[i]; minrank = i;}
    if (globalCountPerPart[i] > max) {max = globalCountPerPart[i]; maxrank = i;}
    sum += globalCountPerPart[i];
  }

  if (me == 0) {
    float avg = (float) sum / (float) checkNparts;
    cout << "Minimum count:  " << min << " on rank " << minrank << endl;
    cout << "Maximum count:  " << max << " on rank " << maxrank << endl;
    cout << "Average count:  " << avg << endl;
    cout << "Total count:    " << sum
         << (sum != origMatrix->getGlobalNumRows()
                 ? "Work was lost; FAIL"
                 : " ")
         << endl;
    cout << "Imbalance:     " << max / avg << endl;
  }

  delete [] countPerPart;
  delete [] globalCountPerPart;

  // Redistribute matrix and vector into new matrix and vector.
  // Can use PartitioningSolution from matrix to redistribute the vectors, too.

  if (me == 0) cout << "Redistributing matrix..." << endl;
  RCP<Matrix_t> redistribMatrix;
  adapter.applyPartitioningSolution(*origMatrix, redistribMatrix,
                                    problem.getSolution());

  if (me == 0) cout << "Redistributing vectors..." << endl;
  RCP<Vector_t> redistribVector;
  MultiVectorAdapter_t adapterVector(origVector); 
  adapterVector.applyPartitioningSolution(*origVector, redistribVector,
                                          problem.getSolution());

  // Create a new product vector for sparse matvec
  RCP<Vector_t> redistribProd;
  redistribProd = Tpetra::createVector<zscalar_t,zlno_t,zgno_t>(
                                       redistribMatrix->getRangeMap());

  // SANITY CHECK
  // check that redistribution is "correct"; perform matvec with
  // original and redistributed matrices/vectors and compare norms.

  if (me == 0) cout << "Matvec original..." << endl;
  origMatrix->apply(*origVector, *origProd);
  zscalar_t origNorm = origProd->norm2();
  if (me == 0)
    cout << "Norm of Original matvec prod:       " << origNorm << endl;

  if (me == 0) cout << "Matvec redistributed..." << endl;
  redistribMatrix->apply(*redistribVector, *redistribProd);
  zscalar_t redistribNorm = redistribProd->norm2();
  if (me == 0)
    cout << "Norm of Redistributed matvec prod:  " << redistribNorm << endl;

  const double epsilon = 0.00000001;
  if (redistribNorm > origNorm+epsilon || redistribNorm < origNorm-epsilon)
    cout << "Mat-Vec product changed; FAIL" << std::endl;
  else
    cout << "PASS" << endl;

  return 0;
}
