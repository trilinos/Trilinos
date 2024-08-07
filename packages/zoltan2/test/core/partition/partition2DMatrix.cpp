// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_MatrixPartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
// #include <Zoltan2_XpetraCrsGraphAdapter.hpp>
// #include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
// #include <iostream>
// #include <limits>
// #include <Teuchos_ParameterList.hpp>
// #include <Teuchos_RCP.hpp>
// #include <Teuchos_FancyOStream.hpp>
// #include <Teuchos_CommandLineProcessor.hpp>
// #include <Tpetra_CrsMatrix.hpp>
// #include <Tpetra_Vector.hpp>
// #include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;

/////////////////////////////////////////////////////////////////////////////
// Program to test Zoltan2 2D partitioning of a TPetra matrix
// (read from a MatrixMarket file) -- modified from partitioning1 test.
// Usage:
//     a.out [--inputFile=filename] [--outputFile=outfile] [--verbose]
// Karen Devine, 2011
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Eventually want to use Teuchos unit tests to vary z2TestLO and
// GO.  For now, we set them at compile time based on whether Tpetra
// is built with explicit instantiation on.  (in Zoltan2_TestHelpers.hpp)

typedef zlno_t z2TestLO;
typedef zgno_t z2TestGO;
typedef zscalar_t z2TestScalar;

typedef Tpetra::CrsMatrix<z2TestScalar, z2TestLO, z2TestGO> SparseMatrix;
typedef Tpetra::CrsGraph<z2TestLO, z2TestGO> SparseGraph;
typedef Tpetra::Vector<z2TestScalar, z2TestLO, z2TestGO> Vector;
typedef Vector::node_type Node;

typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix> SparseMatrixAdapter;
//typedef Zoltan2::TpetraRowMatrixAdapter<SparseMatrix> SparseMatrixAdapter;

//typedef Zoltan2::XpetraCrsGraphAdapter<SparseGraph> SparseGraphAdapter;
typedef Zoltan2::XpetraMultiVectorAdapter<Vector> MultiVectorAdapter;


// Integer vector
typedef Tpetra::Vector<int, z2TestLO, z2TestGO> IntVector;
typedef Zoltan2::XpetraMultiVectorAdapter<IntVector> IntVectorAdapter;

#define epsilon 0.00000001
#define NNZ_IDX 1

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  std::string inputFile = "";        // Matrix Market or Zoltan file to read
  //  std::string outputFile = "";       // Matrix Market or Zoltan file to write
  std::string inputPath = testDataFilePath;  // Directory with input file



  std::string method = "scotch";
  bool verbose = false;              // Verbosity of output
  bool distributeInput = true;
  //  bool haveFailure = false;
  // int nVwgts = 0;
  // int nEwgts = 0;
  int testReturn = 0;

  ////// Establish session.
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("inputPath", &inputPath,
                 "Path to the MatrixMarket or Zoltan file to be read; ",true);

  cmdp.setOption("inputFile", &inputFile,
                 "Name of the Matrix Market or Zoltan file to read; ",true);
  // cmdp.setOption("outputFile", &outputFile,
  //                "Name of the Matrix Market sparse matrix file to write, "
  //                "echoing the input/generated matrix.");
  // cmdp.setOption("method", &method,
  //                "Partitioning method to use:  scotch or parmetis.");
  cmdp.setOption("verbose", "quiet", &verbose,
                 "Print messages and results.");
  // cmdp.setOption("distribute", "no-distribute", &distributeInput,
  //               "indicate whether or not to distribute "
  //               "input across the communicator");


  cmdp.parse(narg, arg);

  RCP<UserInputForTests> uinput;

  // Input file specified; read a matrix
  uinput = rcp(new UserInputForTests(inputPath, inputFile, comm,true, distributeInput));

  RCP<SparseMatrix> origMatrix = uinput->getUITpetraCrsMatrix();

  if (origMatrix->getGlobalNumRows() < 40)
  {
    Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));
    origMatrix->describe(out, Teuchos::VERB_EXTREME);
  }


  if (me == 0)
  {
    std::cout << "NumRows     = " << origMatrix->getGlobalNumRows() << std::endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << std::endl
         << "NumProcs = " << comm->getSize() << std::endl
         << "NumLocalRows (rank 0) = " << origMatrix->getLocalNumRows() << std::endl;
  }


  // ////// Create vectors to use with matrices
  // // Don't use for now, wait until distribution is done.
  // RCP<Vector> origVector, origProd;
  // origProd   = Tpetra::createVector<z2TestScalar,z2TestLO,z2TestGO>(
  //                                   origMatrix->getRangeMap());
  // origVector = Tpetra::createVector<z2TestScalar,z2TestLO,z2TestGO>(
  //                                   origMatrix->getDomainMap());
  // origVector->randomize();

  ////// Specify problem parameters
  Teuchos::ParameterList params;

  params.set("partitioning_approach", "partition");


  // TODO this needs to be added back once parameter is properly validated
  //  params.set("algorithm", "2D Cartesian");


  ////// Create an input adapter for the graph of the Tpetra matrix.
  //  SparseGraphAdapter adapter(origMatrix->getCrsGraph(), nVwgts, nEwgts);

  SparseMatrixAdapter adapter(origMatrix, 0);

  // Zoltan2::TpetraRowMatrixAdapter< User, UserCoord >::TpetraRowMatrixAdapter(const RCP< const User > & inmatrix,
  //                                                                         int nWeightsPerRow = 0
  //                                                                         )

  ////// Create and solve partitioning problem
  Zoltan2::MatrixPartitioningProblem<SparseMatrixAdapter> problem(&adapter, &params);

  try {
    if (me == 0) std::cout << "Calling solve() " << std::endl;

    problem.solve();

    if (me == 0) std::cout << "Done solve() " << std::endl;
  }

  catch (std::runtime_error &e)
  {
    std::cout << "Runtime exception returned from solve(): " << e.what();
    if (!strncmp(e.what(), "BUILD ERROR", 11)) {
      // Catching build errors as exceptions is OK in the tests
      std::cout << " PASS" << std::endl;
      return 0;
    }
    else {
      // All other runtime_errors are failures
      std::cout << " FAIL" << std::endl;
      return -1;
    }
  }
  catch (std::logic_error &e) {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e) {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e) {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }




//   // Don't distribute yet

//   ////// Redistribute matrix and vector into new matrix and vector.
//   if (me == 0) std::cout << "Redistributing matrix..." << std::endl;
//   SparseMatrix *redistribMatrix;
//   SparseMatrixAdapter adapterMatrix(origMatrix);
//   adapterMatrix.applyPartitioningSolution(*origMatrix, redistribMatrix,
//                                           problem.getSolution());
//   if (redistribMatrix->getGlobalNumRows() < 40) {
//     Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));
//     redistribMatrix->describe(out, Teuchos::VERB_EXTREME);
//   }

//   if (me == 0) std::cout << "Redistributing vectors..." << std::endl;
//   Vector *redistribVector;
// //  std::vector<const zscalar_t *> weights;
// //  std::vector<int> weightStrides;
//   MultiVectorAdapter adapterVector(origVector); //, weights, weightStrides);
//   adapterVector.applyPartitioningSolution(*origVector, redistribVector,
//                                           problem.getSolution());

//   RCP<Vector> redistribProd;
//   redistribProd = Tpetra::createVector<z2TestScalar,z2TestLO,z2TestGO>(
//                                        redistribMatrix->getRangeMap());

//   // Test redistributing an integer vector with the same solution.
//   // This test is mostly to make sure compilation always works.
//   RCP<IntVector> origIntVec;
//   IntVector *redistIntVec;
//   origIntVec = Tpetra::createVector<int,z2TestLO,z2TestGO>(
//                                         origMatrix->getRangeMap());
//   for (size_t i = 0; i < origIntVec->getLocalLength(); i++)
//     origIntVec->replaceLocalValue(i, me);

//   IntVectorAdapter int_vec_adapter(origIntVec);
//   int_vec_adapter.applyPartitioningSolution(*origIntVec, redistIntVec,
//                                              problem.getSolution());
//   int origIntNorm = origIntVec->norm1();
//   int redistIntNorm = redistIntVec->norm1();
//   if (me == 0) std::cout << "IntegerVectorTest:  " << origIntNorm << " == "
//                          << redistIntNorm << " ?";
//   if (origIntNorm != redistIntNorm) {
//     if (me == 0) std::cout << " FAIL" << std::endl;
//     haveFailure = true;
//   }
//   else if (me == 0) std::cout << " OK" << std::endl;
//   delete redistIntVec;

//   ////// Verify that redistribution is "correct"; perform matvec with
//   ////// original and redistributed matrices/vectors and compare norms.

//   if (me == 0) std::cout << "Matvec original..." << std::endl;
//   origMatrix->apply(*origVector, *origProd);
//   z2TestScalar origNorm = origProd->norm2();
//   if (me == 0)
//     std::cout << "Norm of Original matvec prod:       " << origNorm << std::endl;

//   if (me == 0) std::cout << "Matvec redistributed..." << std::endl;
//   redistribMatrix->apply(*redistribVector, *redistribProd);
//   z2TestScalar redistribNorm = redistribProd->norm2();
//   if (me == 0)
//     std::cout << "Norm of Redistributed matvec prod:  " << redistribNorm << std::endl;

//   if (redistribNorm > origNorm+epsilon || redistribNorm < origNorm-epsilon) {
//     testReturn = 1;
//     haveFailure = true;
//   }

//   delete redistribVector;
//   delete redistribMatrix;

//   if (me == 0) {
//     if (testReturn) {
//       std::cout << "Mat-Vec product changed; FAIL" << std::endl;
//       haveFailure = true;
//     }
//     if (!haveFailure)
//       std::cout << "PASS" << std::endl;
//   }

  std::cout << "Finished" << std::endl;

  return testReturn;
}
