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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::OSTab;
using Teuchos::ParameterList;
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

const bool callFillComplete = true;
const bool tolerant = false;

namespace {

const char matrix_symRealSmall[] = 
"%%MatrixMarket matrix coordinate real general\n"
"5 5 13\n"
"1 1  2.0000000000000e+00\n"
"1 2  -1.0000000000000e+00\n"
"2 1  -1.0000000000000e+00\n"
"2 2  2.0000000000000e+00\n"
"2 3  -1.0000000000000e+00\n"
"3 2  -1.0000000000000e+00\n"
"3 3  2.0000000000000e+00\n"
"3 4  -1.0000000000000e+00\n"
"4 3  -1.0000000000000e+00\n"
"4 4  2.0000000000000e+00\n"
"4 5  -1.0000000000000e+00\n"
"5 4  -1.0000000000000e+00\n"
"5 5  2.0000000000000e+00\n";

} // namespace (anonymous)


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase0, ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  // Whether to print copious debugging output to stderr when doing
  // Matrix Market input and output. 
  const bool debug = false;

  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;

  out << "Test: CrsMatrix Matrix Market I/O, w/ Map with index base 0" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  // Forestall compiler warnings about unused variables.
  (void) myRank;
  (void) numProcs;

  out << "Original sparse matrix:" << endl;
  out << matrix_symRealSmall << endl;

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  const GO indexBase = 0;
  RCP<const map_type> rowMap = 
    rcp (new map_type (globalNumElts, indexBase, comm, 
		       Tpetra::GloballyDistributed, node));

  out << "Reading in the matrix" << endl;
  std::istringstream inStr (matrix_symRealSmall);
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<crs_matrix_type> A = 
    reader_type::readSparse (inStr, rowMap, colMap, domainMap, rangeMap, 
			     callFillComplete, tolerant, debug);

  out << "Writing out the matrix" << endl;
  std::ostringstream outStr;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparse (out, A, debug);
  
  out << "Result of writing the CrsMatrix:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << outStr.str () << endl;
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase1, ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  // Whether to print copious debugging output to stderr when doing
  // Matrix Market input and output. 
  const bool debug = true;

  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;

  out << "Test: CrsMatrix Matrix Market I/O, w/ Map with index base 1" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  // Forestall compiler warnings about unused variables.
  (void) myRank;
  (void) numProcs;

  out << "Original sparse matrix:" << endl;
  out << matrix_symRealSmall << endl;

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  const GO indexBase = 1;
  RCP<const map_type> rowMap = 
    rcp (new map_type (globalNumElts, indexBase, comm, 
		       Tpetra::GloballyDistributed, node));

  out << "Reading in the matrix" << endl;
  std::istringstream inStr (matrix_symRealSmall);
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<crs_matrix_type> A = 
    reader_type::readSparse (inStr, rowMap, colMap, domainMap, rangeMap, 
			     callFillComplete, tolerant, debug);

  out << "Writing out the matrix" << endl;
  std::ostringstream outStr;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparse (out, A, debug);
  
  out << "Result of writing the CrsMatrix:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << outStr.str () << endl;
  }
}


// Unit test macro isn't smart enough to deal with namespace qualifications.
typedef Kokkos::DefaultNode::DefaultNodeType the_node_type;

// We instantiate tests for all combinations of the following parameters:
// - indexBase = {0, 1}
// - ST = {double, float}
// - GO = {int, long}
//
// We should really use the Tpetra ETI system to control which GO we
// test here, but int and long are the two most important cases.

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, double, int, int, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, double, int, int, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, double, int, long, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, double, int, long, the_node_type )

