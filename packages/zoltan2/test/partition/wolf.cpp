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
#include <Zoltan2_PartitioningProblem.hpp>
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
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;
using namespace std;

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

//typedef Tpetra::MultiVector<z2TestScalar, z2TestLO, z2TestGO, znode_t> tMVector_t;
typedef Tpetra::MultiVector<z2TestScalar, z2TestLO, z2TestGO,znode_t> tMVector_t;


//typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix_t> SparseMatrixAdapter;
typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix_t,tMVector_t> SparseMatrixAdapter_t;

typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> MultiVectorAdapter_t;




#define epsilon 0.00000001

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string inputFile = "";        // Matrix Market or Zoltan file to read
  std::string inputPath = testDataFilePath;  // Directory with input file
  bool distributeInput = true;
  int success = 0;
  int numParts = 2;


  //////////////////////////////////////////////////////////////////////
  ////// Establish session.
  //////////////////////////////////////////////////////////////////////
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  int me = comm->getRank();
  //////////////////////////////////////////////////////////////////////

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
    parseReturn= cmdp.parse( argc, argv );

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
    cout << "NumRows     = " << origMatrix->getGlobalNumRows() << endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << endl
         << "NumProcs = " << comm->getSize() << endl;
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

  //MMW add weights later

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("partitioning_approach", "partition");
  params.set("algorithm", "wolf");
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
  Zoltan2::PartitioningProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);

  try
  {
    if (me == 0) cout << "Calling solve() " << endl;
    problem.solve();
    if (me == 0) cout << "Done solve() " << endl;
  }
  catch (std::runtime_error &e)
  {
    cout << "Runtime exception returned from solve(): " << e.what();
    if (!strncmp(e.what(), "BUILD ERROR", 11)) {
      // Catching build errors as exceptions is OK in the tests
      cout << " PASS" << endl;
      return 0;
    }
    else {
      // All other runtime_errors are failures
      cout << " FAIL" << endl;
      return -1;
    }
  }
  catch (std::logic_error &e)
  {
    cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << endl;
    return -1;
  }
  catch (std::exception &e)
  {
    cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << endl;
    return -1;
  }

  //////////////////////////////////////////////////////////////////////

  delete ca;


  std::cout << "PASS" << std::endl;
  return success;
}
////////////////////////////////////////////////////////////////////////////////

