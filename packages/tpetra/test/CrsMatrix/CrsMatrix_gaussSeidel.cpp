/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_ETIHelperMacros.h>

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

//
// Tests for Tpetra::CrsMatrix::gaussSeidel().
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, gaussSeidelSerial, LocalOrdinalType, GlobalOrdinalType, ScalarType, NodeType )
{
  using Tpetra::createContigMapWithNode;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createMultiVector;
  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_const_cast;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::OrdinalTraits;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_SUM;
  using Teuchos::REDUCE_MIN;
  using Teuchos::ScalarTraits;
  using Teuchos::tuple;
  using Teuchos::TypeNameTraits;
  using std::cerr;
  using std::endl;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef NodeType node_type;

  // Typedefs derived from the above canonical typedefs.
  typedef ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef ScalarTraits<magnitude_type> STM;
  typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // Abbreviation typedefs.
  typedef scalar_type ST;
  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
  typedef node_type NT;

  // CrsMatrix specialization to use in this test.
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;

  // CrsGraph specialization corresponding to crs_matrix_type (the
  // CrsMatrix specialization).
  typedef Tpetra::CrsGraph<LO, GO, NT, typename crs_matrix_type::mat_solve_type> crs_graph_type;

  // MultiVector specialization corresponding to crs_matrix_type.
  typedef Tpetra::MultiVector<ST, LO, GO, NT> multivector_type;
  // Vector specialization corresponding to crs_matrix_type.
  typedef Tpetra::Vector<ST, LO, GO, NT> vector_type;


  ////////////////////////////////////////////////////////////////////  
  // HERE BEGINS THE TEST.
  ////////////////////////////////////////////////////////////////////

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

  // Get the default communicator.
  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();

  if (myRank == 0) {
    out << "Test with " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;
  }

#if 0
  // This test doesn't make much sense if there is only one MPI
  // process.  We let it pass trivially in that case.
  if (numProcs == 1) {
    out << "Number of processes in world is one; test passes trivially." << endl;
    return;
  }
#endif // 0

  // Get a Kokkos Node instance.  It would be nice if we could pass in
  // parameters here, but threads don't matter for this test, since
  // the current Gauss-Seidel implementation doesn't use them.

  if (myRank == 0) {
    out << "Creating Kokkos Node of type " << TypeNameTraits<node_type>::name () << endl;
  }
  RCP<node_type> node;
  {
    ParameterList pl; // Kokkos Node types require a PL inout.
    node = rcp (new node_type (pl));
  }

  // Let's build ourselves the graph of a 2-D Laplacian (Dirichlet
  // boundary conditions) on a square numGlobalPoints x
  // numGlobalPoints mesh.  Each process gets a numLocalPoints x
  // numGlobalPoints slab.

  const LO numLocalPoints = 9;
  const GO numGlobalPoints = numProcs * numLocalPoints;
  
  // Number of rows in the matrix owned by each process.
  const LO numLocalRows = numLocalPoints * numGlobalPoints;

  // Number of (global) rows and columns in the matrix.
  const GO numGlobalRows = numGlobalPoints * numGlobalPoints;
  const GO numGlobalCols = numGlobalRows;
  // Prevent compile warning for unused variable.
  // (It's not really "variable" if it's const, but oh well.)
  (void) numGlobalCols; 

  if (myRank == 0) {
    out << "Creating contiguous row Map" << endl;
  }

  // Create a contiguous row Map, with numLocalRows rows per process.
  RCP<const map_type> rowMap = createContigMapWithNode<LO, GO, NT> (INVALID, numLocalRows, comm, node);

  // The Gauss-Seidel kernel requires that the row, domain, and range
  // Maps all be the same.
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;

  // Min and max row index of this process.
  const GO globalMinRow = rowMap->getMinGlobalIndex ();
  const GO globalMaxRow = rowMap->getMaxGlobalIndex ();
  const GO globalMinAllRow = rowMap->getMinAllGlobalIndex ();
  const GO globalMaxAllRow = rowMap->getMaxAllGlobalIndex ();

  if (myRank == 0) {
    out << "Creating graph" << endl;
  }

  // Create a numGlobalRows by numGlobalCols graph and set its
  // structure.  Every process sets its diagonal entries (which it
  // owns), and entries 1 and numGlobalPoints away in each direction
  // (if valid).
  RCP<const crs_graph_type> graph;
  {
    // We have a good upper bound for the number of entries per row, so use static profile.
    RCP<crs_graph_type> nonconstGraph (new crs_graph_type (rowMap, 5, Tpetra::StaticProfile));

    for (GO globalRow = globalMinRow; globalRow <= globalMaxRow; ++globalRow) {
      Teuchos::Array<GO> indices;
      if (globalRow - numGlobalPoints >= globalMinAllRow) {
	indices.push_back (globalRow - numGlobalPoints);
      }
      if (globalRow - 1 >= globalMinAllRow) {
	indices.push_back (globalRow - 1);
      }
      indices.push_back (globalRow);
      if (globalRow + 1 <= globalMaxAllRow) {
	indices.push_back (globalRow + 1);
      }
      if (globalRow + numGlobalPoints <= globalMaxAllRow) {
	indices.push_back (globalRow + numGlobalPoints);
      }
      nonconstGraph->insertGlobalIndices (globalRow, indices ());
    }

    nonconstGraph->fillComplete (domainMap, rangeMap);
    graph = rcp_const_cast<const crs_graph_type> (nonconstGraph);
  }

  if (myRank == 0) {
    out << "Creating matrix" << endl;
  }

  // Create the matrix, using the above graph.
  RCP<crs_matrix_type> matrix (new crs_matrix_type (graph));
  {
    for (GO globalRow = globalMinRow; globalRow <= globalMaxRow; ++globalRow) {
      Teuchos::Array<GO> indices;
      Teuchos::Array<ST> values;
      if (globalRow - numGlobalPoints >= globalMinAllRow) {
	indices.push_back (globalRow - numGlobalPoints);
	values.push_back (-STS::one ());
      }
      if (globalRow - 1 >= globalMinAllRow) {
	indices.push_back (globalRow - 1);
	values.push_back (-STS::one ());
      }
      indices.push_back (globalRow);
      values.push_back (as<ST> (4));
      if (globalRow + 1 <= globalMaxAllRow) {
	indices.push_back (globalRow + 1);
	values.push_back (-STS::one ());
      }
      if (globalRow + numGlobalPoints <= globalMaxAllRow) {
	indices.push_back (globalRow + numGlobalPoints);
	values.push_back (-STS::one ());
      }
      matrix->replaceGlobalValues (globalRow, indices (), values ());
    }

    if (myRank == 0) {
      out << "Calling fillComplete on the matrix" << endl;
    }
    matrix->fillComplete (domainMap, rangeMap);
  }

  const bool dumpMatrix = false;
  if (dumpMatrix) {
    const std::string filename ("A.mtx");
    const std::string scalarType = Teuchos::TypeNameTraits<ST>::name ();
    const std::string matName = "A-" + scalarType;
    typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
    writer_type::writeSparseFile (filename, matrix, matName, "Gauss-Seidel test matrix");
  }

  if (myRank == 0) {
    cerr << "Extracting diagonal" << endl;
  }
  RCP<vector_type> D = rcp (new vector_type (rowMap));
  matrix->getLocalDiagCopy (*D);

  if (myRank == 0) {
    cerr << "Making vectors" << endl;
  }

  // Make (multi)vectors for the initial guess (which will also be the
  // solution vector), the right-hand side, and the residual.  We're
  // only testing a single right-hand side for now.  Make all the
  // vectors first in the column Map, with the actual vector given to
  // Gauss-Seidel in the domain / range Map.

  RCP<const map_type> colMap = graph->getColMap ();
  RCP<multivector_type> X_colMap = createMultiVector<ST, LO, GO, NT> (colMap, 1);
  RCP<multivector_type> X_exact_colMap = createMultiVector<ST, LO, GO, NT> (colMap, 1);
  RCP<multivector_type> B_colMap = createMultiVector<ST, LO, GO, NT> (colMap, 1);
  RCP<multivector_type> R_colMap = createMultiVector<ST, LO, GO, NT> (colMap, 1);
  RCP<multivector_type> X = X_colMap->offsetViewNonConst (domainMap, 0);
  RCP<multivector_type> X_exact = X_exact_colMap->offsetViewNonConst (domainMap, 0);
  RCP<multivector_type> B = B_colMap->offsetViewNonConst (rangeMap, 0);
  RCP<multivector_type> R = R_colMap->offsetViewNonConst (rangeMap, 0);

  // Set the exact solution and right-hand side.
  X_exact->randomize ();
  matrix->apply (*X_exact, *B);

  const int maxNumIters = 10;
  Array<magnitude_type> residNorms (maxNumIters);

  if (myRank == 0) {
    cerr << "Iterating" << endl;
  }

  int localSuccess = 1;
  int globalSuccess = 1;
  int smallestFailingRank = 0;
  for (int iter = 0; iter < maxNumIters; ++iter) {
    // R = B - A * X
    matrix->apply (*X, *R); // R = A * X
    R->update (STS::one(), *B, -STS::one()); // R = 1*B - 1*R
    // residNorms[iter] = \|R\|_2
    Array<magnitude_type> normTemp (1);
    R->norm2 (normTemp);
    residNorms[iter] = normTemp[0];

    std::string exMsg;
    try {
      matrix->gaussSeidel (*B, *X, *D, STS::one(), Tpetra::Symmetric, 1);
    }
    catch (std::exception& e) {
      exMsg = e.what ();
      localSuccess = 0;
    }
    reduceAll (*comm, REDUCE_SUM, localSuccess, outArg (globalSuccess));
    reduceAll (*comm, REDUCE_SUM, localSuccess, outArg (globalSuccess));
    if (globalSuccess < numProcs) {
      // Compute min rank that failed.  For procs that didn't fail,
      // set their "rank" to numProcs.
      const int inRank = (localSuccess == 1) ? numProcs : myRank;
      reduceAll (*comm, REDUCE_MIN, inRank, outArg (smallestFailingRank));
      // The min-failing-rank process prints the error message.
      if (myRank == smallestFailingRank) {
	cerr << "Proc " << myRank << " (possibly among others) failed.  "
	     << "Here is its exception message: " << exMsg;
      }
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "");
    }      
    if (myRank == 0) {
      cerr << "Iteration " << iter+1 << " of " << maxNumIters 
	   << ": ||R||_2 = " << normTemp[0] << endl;
    }
  }

#if 0
  bool globalSuccess = true;
  {
    int globalSuccess_int = 1;
    reduceAll (*comm, Teuchos::REDUCE_MIN, localSuccess ? 1 : 0, outArg (globalSuccess_int));
    globalSuccess = (globalSuccess_int != 0);
  }

  if (! globalSuccess) {
    // Print out the failure messages on all processes.
    for (int p = 0; p < numProcs; ++p) {
      if (p == myRank) {
	out << failMsg.str () << endl;
	out << "Proc " << myRank << ": localSuccess = " << localSuccess << ", globalSuccess = " << globalSuccess << endl;
	//      std::flush (out);
      }
      // Do some barriers to allow output to finish.
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    }
  }

  TEST_EQUALITY_CONST(globalSuccess, true);
#endif // 0
}

//////////////////////////////////////////////////////////////////////
// INSTANTIATE THE TEMPLATED UNIT TESTS
//////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, gaussSeidelSerial, LO, GO, SCALAR, NODE )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )



