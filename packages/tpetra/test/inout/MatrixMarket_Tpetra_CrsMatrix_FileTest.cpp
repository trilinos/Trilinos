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
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Util.hpp> // sort2
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

namespace {

std::string mapOutputFilename;
std::string matrixOutputFilename;
std::string mapInputFilename;
std::string matrixInputFilename;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
  clp.setOption ("mapOutputFilename", &mapOutputFilename,
                 "Name of file to which to write temporary Map output");
  clp.setOption ("matrixOutputFilename", &matrixOutputFilename,
                 "Name of file to which to write temporary CrsMatrix output");
  clp.setOption ("mapInputFilename", &mapInputFilename,
                 "Name of file from which to read the initial Map");
  clp.setOption ("matrixInputFilename", &matrixInputFilename,
                 "Name of file from which to read the initial CrsMatrix "
                 "(using the input Map as its row Map)");
}

// Input matrices must be fill complete.
template<class CrsMatrixType>
bool
compareCrsMatrix (const CrsMatrixType& A_orig, const CrsMatrixType& A)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;

  if (! A_orig.getRowMap ()->isSameAs (* (A.getRowMap ()))) {
    return false;
  }
  else if (! A_orig.getColMap ()->isSameAs (* (A.getColMap ()))) {
    return false;
  }
  else if (! A_orig.getDomainMap ()->isSameAs (* (A.getDomainMap ()))) {
    return false;
  }
  else if (! A_orig.getRangeMap ()->isSameAs (* (A.getRangeMap ()))) {
    return false;
  }
  else {
    //
    // Are my local matrices equal?
    //
    RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
    int localEqual = 1;

    Array<GO> indOrig, ind;
    Array<ST> valOrig, val;
    size_t numEntriesOrig = 0;
    size_t numEntries = 0;

    ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
    const size_type numLocalElts = localElts.size ();
    for (size_type k = 0; k < numLocalElts; ++k) {
      const GO globalRow = localElts[k];
      numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
      numEntries = A.getNumEntriesInGlobalRow (globalRow);

      if (numEntriesOrig != numEntries) {
        localEqual = 0;
        break;
      }
      indOrig.resize (numEntriesOrig);
      valOrig.resize (numEntriesOrig);
      A_orig.getGlobalRowCopy (globalRow, indOrig (), valOrig (), numEntriesOrig);
      ind.resize (numEntries);
      val.resize (numEntries);
      A.getGlobalRowCopy (globalRow, ind (), val (), numEntries);

      // Global row entries are not necessarily sorted.  Sort them so
      // we can compare them.
      Tpetra::sort2 (indOrig.begin (), indOrig.end (), valOrig.begin ());
      Tpetra::sort2 (ind.begin (), ind.end (), val.begin ());

      for (size_t entryIndex = 0; entryIndex < numEntries; ++entryIndex) {
        // Values should be _exactly_ equal.
        if (indOrig[k] != ind[k] || valOrig[k] != val[k]) {
          localEqual = 0;
          break;
        }
      }
    }

    int globalEqual = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
    return globalEqual == 1;
  }
}


void
testReadAndWriteFile (Teuchos::FancyOStream& out,
                      const std::string& mapOutFile,
                      const std::string& matrixOutFile,
                      const std::string& mapInFile,
                      const std::string& matrixInFile)
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  typedef double ST;
  typedef int LO;
  typedef long GO;
  typedef Kokkos::SerialNode NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  const bool tolerant = false;
  // Whether to print copious debugging output to stderr when doing
  // Matrix Market input and output.
  const bool debug = false;

  bool result = true; // current Boolean result; reused below
  bool success = true; // used by TEST_EQUALITY

  out << "Test: Matrix Market, CrsMatrix from row Map" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

  out << "Reading the (row) Map" << endl;
  RCP<const map_type> rowMap =
    reader_type::readMapFile (mapInFile, comm, node, tolerant, debug);

  out << "Writing the Map" << endl;
  writer_type::writeMapFile (mapOutFile, *rowMap);

  out << "Reading the written Map back in" << endl;
  RCP<const map_type> rowMap_out =
    reader_type::readMapFile (mapOutFile, comm, node, tolerant, debug);

  out << "Comparing the two Maps" << endl;
  result = rowMap->isSameAs (*rowMap_out);
  TEST_EQUALITY( result, true );

  out << "Reading the CrsMatrix" << endl;
  const bool callFillComplete = true;
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  RCP<crs_matrix_type> A_in =
    reader_type::readSparseFile (matrixInFile,
                                 rowMap, colMap, domainMap, rangeMap,
                                 callFillComplete, tolerant, debug);

  out << "Writing the CrsMatrix" << endl;
  writer_type::writeSparseFile (matrixOutFile, A_in, debug);

  out << "Reading the written CrsMatrix back in" << endl;
  RCP<crs_matrix_type> A_out =
    reader_type::readSparseFile (matrixOutFile,
                                 rowMap, colMap, domainMap, rangeMap,
                                 callFillComplete, tolerant, debug);

  out << "Test the two matrices for equality" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_in, *A_out);
  TEST_EQUALITY( result, true );
}

} // namespace (anonymous)

TEUCHOS_UNIT_TEST( MatrixMarket, CrsMatrixFileTest )
{
  testReadAndWriteFile (out, mapOutputFilename, matrixOutputFilename,
                        mapInputFilename, matrixInputFilename);
}

