// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <iostream>
#include <limits>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;

/////////////////////////////////////////////////////////////////////////////
// Eventually want to use Teuchos unit tests to vary z2TestLO and
// GO.  For now, we set them at compile time based on whether Tpetra
// is built with explicit instantiation on.  (in Zoltan2_TestHelpers.hpp)

typedef zlno_t z2TestLO;
typedef zgno_t z2TestGO;
typedef zscalar_t z2TestScalar;

typedef Tpetra::CrsMatrix<z2TestScalar, z2TestLO, z2TestGO> SparseMatrix_t;
typedef Tpetra::Vector<z2TestScalar, z2TestLO, z2TestGO> Vector;
typedef Vector::node_type Node;

typedef Tpetra::MultiVector<z2TestScalar, z2TestLO, z2TestGO,znode_t> tMVector_t;


//typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix_t> SparseMatrixAdapter;
typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix_t,tMVector_t> SparseMatrixAdapter_t;

typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> MultiVectorAdapter_t;

#define epsilon 0.00000001

int testNDwithRCB(RCP<SparseMatrix_t> &origMatrix,RCP<tMVector_t> &coords, int numParts, int me);
int testNDwithPHG(RCP<SparseMatrix_t> &origMatrix,int numParts, int me);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  //////////////////////////////////////////////////////////////////////
  ////// Establish session.
  //////////////////////////////////////////////////////////////////////
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  //////////////////////////////////////////////////////////////////////

  std::string inputFile = "";        // Matrix Market or Zoltan file to read
  std::string inputPath = testDataFilePath;  // Directory with input file
  bool distributeInput = true;
  int success = 0;
  int numParts = 2;

  //////////////////////////////////////////////////////////////////////
  // Read run-time options.
  //////////////////////////////////////////////////////////////////////
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("inputPath", &inputPath,
                 "Path to the MatrixMarket or Zoltan file to be read; "
                 "if not specified, a default path will be used.");
  cmdp.setOption("inputFile", &inputFile,
                 "Name of the Matrix Market or Zoltan file to read; "
                 "");
  cmdp.setOption("distribute", "no-distribute", &distributeInput,
                "for Zoltan input files only, "
                "indicate whether or not to distribute "
                "input across the communicator");
  cmdp.setOption("numParts", &numParts,
                 "Global number of parts;");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn= cmdp.parse( narg, arg );

  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
  {
    return 0;
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Construct matrix from file
  //////////////////////////////////////////////////////////////////////
  RCP<UserInputForTests> uinput;

  if (inputFile != "")   // Input file specified; read a matrix
  {
    uinput = rcp(new UserInputForTests(inputPath, inputFile, comm,
                                       true, distributeInput));
  }
  else
  {
    std::cout << "Input file must be specified." << std::endl;
  }

  RCP<SparseMatrix_t> origMatrix = uinput->getUITpetraCrsMatrix();


  if (me == 0)
  {
    std::cout << "NumRows     = " << origMatrix->getGlobalNumRows() << std::endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << std::endl
         << "NumProcs = " << comm->getSize() << std::endl
         << "NumParts = " << numParts << std::endl;
  }

  if (origMatrix->getGlobalNumRows() < 40)
  {
    Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));
    origMatrix->describe(out, Teuchos::VERB_EXTREME);
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Get coordinates, corresponding to graph vertices
  //////////////////////////////////////////////////////////////////////
  RCP<tMVector_t> coords;
  try
  {
    coords = uinput->getUICoordinates();
  }
  catch(...)
  {
    if (me == 0)
      std::cout << "FAIL: get coordinates" << std::endl;
    return 1;
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Test ND ordering with RCB to compute the separator
  //////////////////////////////////////////////////////////////////////
  testNDwithRCB(origMatrix,coords,numParts,me);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Test ND ordering with PHG to compute the separator
  //////////////////////////////////////////////////////////////////////
  testNDwithPHG(origMatrix,numParts,me);
  //////////////////////////////////////////////////////////////////////



  return success;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int testNDwithRCB(RCP<SparseMatrix_t> &origMatrix,RCP<tMVector_t> &coords, int numParts, int me)
{
  int success=0;

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("order_method", "nd");
  params.set("edge_separator_method", "rcb");
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create an input adapter for the Tpetra matrix.
  //////////////////////////////////////////////////////////////////////
  SparseMatrixAdapter_t matAdapter(origMatrix);

  MultiVectorAdapter_t *ca = NULL;

  try
  {
       ca = new MultiVectorAdapter_t(coords);
  }
  catch(...)
  {
    if (me == 0)
      std::cout << "FAIL: vector adapter" << std::endl;
    return 1;
  }

  matAdapter.setCoordinateInput(ca);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create and solve partitioning problem
  //////////////////////////////////////////////////////////////////////
  Zoltan2::OrderingProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);


  try
  {
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
  catch (std::logic_error &e)
  {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e)
  {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }

  //////////////////////////////////////////////////////////////////////

  delete ca;


  std::cout << "PASS" << std::endl;
  return success;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int testNDwithPHG(RCP<SparseMatrix_t> &origMatrix,int numParts, int me)
{
  int success=0;

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("order_method", "nd");
  params.set("edge_separator_method", "phg");
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create an input adapter for the Tpetra matrix.
  //////////////////////////////////////////////////////////////////////
  SparseMatrixAdapter_t matAdapter(origMatrix);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create and solve partitioning problem
  //////////////////////////////////////////////////////////////////////
  Zoltan2::OrderingProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);


  try
  {
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
  catch (std::logic_error &e)
  {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e)
  {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  //////////////////////////////////////////////////////////////////////

  std::cout << "PASS" << std::endl;
  return success;
}
////////////////////////////////////////////////////////////////////////////////
