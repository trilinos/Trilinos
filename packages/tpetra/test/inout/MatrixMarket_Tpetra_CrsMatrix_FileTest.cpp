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
compareCrsMatrixMaps (const CrsMatrixType& A_orig, const CrsMatrixType& A, Teuchos::FancyOStream& out)
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

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));

  bool globalAllSame = true;
  if (! A_orig.getRowMap ()->isSameAs (* (A.getRowMap ()))) {
    out << "Row Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getColMap ()->isSameAs (* (A.getColMap ()))) {
    out << "Column Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getDomainMap ()->isSameAs (* (A.getDomainMap ()))) {
    out << "Domain Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getRangeMap ()->isSameAs (* (A.getRangeMap ()))) {
    out << "Range Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (globalAllSame) {
    out << "All Maps are the same" << endl;
  }
  return globalAllSame;
}

// Input matrices must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsMatrixType>
bool
compareCrsMatrix (const CrsMatrixType& A_orig, const CrsMatrixType& A, Teuchos::FancyOStream& out)
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

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
  int localEqual = 1;

  //
  // Are my local matrices equal?
  //
  Array<GO> indOrig, ind;
  Array<ST> valOrig, val;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
  const size_type numLocalElts = localElts.size ();
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
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

    for (size_t k = 0; k < numEntries; ++k) {
      // Values should be _exactly_ equal.
      if (indOrig[k] != ind[k] || valOrig[k] != val[k]) {
        localEqual = 0;
        break;
      }
    }
  }

  RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
  int globalEqual = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
  return globalEqual == 1;
}

template<class IT1, class IT2>
void
merge2 (IT1& indResultOut, IT2& valResultOut, IT1 indBeg, IT1 indEnd, IT2 valBeg, IT2 valEnd)
{
  if (indBeg == indEnd) {
    indResultOut = indBeg; // It's allowed for indResultOut to alias indEnd
    valResultOut = valBeg; // It's allowed for valResultOut to alias valEnd
  }
  else {
    IT1 indResult = indBeg;
    IT2 valResult = valBeg;
    if (indBeg != indEnd) {
      ++indBeg;
      ++valBeg;
      while (indBeg != indEnd) {
        if (*indResult == *indBeg) { // adjacent column indices equal
          *valResult += *valBeg; // merge entries by adding their values together
        } else { // adjacent column indices not equal
          *(++indResult) = *indBeg; // shift over the index
          *(++valResult) = *valBeg; // shift over the value
        }
        ++indBeg;
        ++valBeg;
      }
      ++indResult; // exclusive end of merged result
      ++valResult; // exclusive end of merged result
      indEnd = indResult;
      valEnd = valResult;
    }
    indResultOut = indResult;
    valResultOut = valResult;
  }
}

template<class IT1, class IT2, class BinaryFunction>
void
merge2 (IT1& indResultOut, IT2& valResultOut,
        IT1 indBeg, IT1 indEnd,
        IT2 valBeg, IT2 valEnd,
        BinaryFunction f)
{
  if (indBeg == indEnd) {
    indResultOut = indBeg; // It's allowed for indResultOut to alias indEnd
    valResultOut = valBeg; // It's allowed for valResultOut to alias valEnd
  }
  else {
    IT1 indResult = indBeg;
    IT2 valResult = valBeg;
    if (indBeg != indEnd) {
      ++indBeg;
      ++valBeg;
      while (indBeg != indEnd) {
        if (*indResult == *indBeg) { // adjacent column indices equal
          *valResult = f (*valResult, *valBeg); // merge entries by adding their values together
        } else { // adjacent column indices not equal
          *(++indResult) = *indBeg; // shift over the index
          *(++valResult) = *valBeg; // shift over the value
        }
        ++indBeg;
        ++valBeg;
      }
      ++indResult; // exclusive end of merged result
      ++valResult; // exclusive end of merged result
      indEnd = indResult;
      valEnd = valResult;
    }
    indResultOut = indResult;
    valResultOut = valResult;
  }
}

// Input matrices must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsMatrixType>
bool
compareCrsMatrixValues (const CrsMatrixType& A_orig,
                        const CrsMatrixType& A,
                        Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  using std::endl;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));

  //
  // Are my local matrices equal?
  //
  Array<GO> indOrig, ind;
  Array<ST> valOrig, val;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
  const size_type numLocalElts = localElts.size ();
  MT localDiff = STM::zero (); // \sum_{i,j} |A_orig(i,j) - A(i,j)| locally
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    indOrig.resize (numEntriesOrig);
    valOrig.resize (numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig (), valOrig (), numEntriesOrig);
    ind.resize (numEntries);
    val.resize (numEntries);
    A.getGlobalRowCopy (globalRow, ind (), val (), numEntries);

    // Global row entries are not necessarily sorted.  Sort them
    // (and their values with them) so we can merge their values.
    Tpetra::sort2 (indOrig.begin (), indOrig.end (), valOrig.begin ());
    Tpetra::sort2 (ind.begin (), ind.end (), val.begin ());

    //
    // Merge repeated values in each set of indices and values.
    //

    typename Array<GO>::iterator indOrigIter = indOrig.begin ();
    typename Array<ST>::iterator valOrigIter = valOrig.begin ();
    typename Array<GO>::iterator indOrigEnd = indOrig.end ();
    typename Array<ST>::iterator valOrigEnd = valOrig.end ();
    merge2 (indOrigEnd, valOrigEnd, indOrigIter, indOrigEnd, valOrigIter, valOrigEnd);

    typename Array<GO>::iterator indIter = ind.begin ();
    typename Array<ST>::iterator valIter = val.begin ();
    typename Array<GO>::iterator indEnd = ind.end ();
    typename Array<ST>::iterator valEnd = val.end ();
    merge2 (indEnd, valEnd, indIter, indEnd, valIter, valEnd);

    //
    // Compare the merged sets of entries.
    //

    indOrigIter = indOrig.begin ();
    indIter = ind.begin ();
    valOrigIter = valOrig.begin ();
    valIter = val.begin ();
    while (indOrigIter != indOrigEnd && indIter != indEnd) {
      const GO j_orig = *indOrigIter;
      const GO j = *indIter;

      if (j_orig < j) { // entry is in A_orig, not in A
        localDiff += STS::magnitude (*valOrigIter);
        ++indOrigIter;
        ++valOrigIter;
      } else if (j_orig > j) { // entry is in A, not A_orig
        localDiff += STS::magnitude (*valIter);
        ++indIter;
        ++valIter;
      } else { // j_orig == j: entry is in both matrices
        localDiff += STS::magnitude (*valOrigIter - *valIter);
        ++indOrigIter;
        ++valOrigIter;
        ++indIter;
        ++valIter;
      }
    }
  }

  RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
  const int myRank = comm->getRank ();
  std::ostringstream os;
  os << "Values are ";
  if (localDiff == STM::zero ()) {
    os << "the same on process " << myRank << endl;
  } else {
    os << "NOT the same on process " << myRank
       << ": \\sum_{i,j} |A_orig(i,j) - A(i,j)| = " << localDiff << endl;
  }
  out << os.str ();

  int globalEqual = 0;
  int localEqual = (localDiff == STM::zero ()) ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
  return globalEqual == 1;
}


bool
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

  out << "Test the two matrices for exact equality of structure and values" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_in, *A_out, out);
  {
    OSTab tab (rcpFromRef (out));
    out << "- Matrices are " << (result ? "" : "NOT ") << "equal" << endl;
  }
  TEST_EQUALITY( result, true );

  out << "Test the two matrices for exact equality of values, independent of local structure" << endl;
  result = compareCrsMatrixValues<crs_matrix_type> (*A_in, *A_out, out);
  TEST_EQUALITY( result, true );

  return result;
}

} // namespace (anonymous)


TEUCHOS_UNIT_TEST( MatrixMarket, Merge2 )
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::toString;
  using std::endl;
  typedef Array<int>::size_type size_type;
  const size_type origNumEntries = 8;

  Array<int> ind (origNumEntries);
  ind[0] =  0;
  ind[1] =  1;
  ind[2] =  1;
  ind[3] =  3;
  ind[4] = -1;
  ind[5] = -1;
  ind[6] = -1;
  ind[7] =  0;
  Array<int> indCopy = ind; // deep copy

  Array<double> val (origNumEntries);
  val[0] =  42.0;
  val[1] =  -4.0;
  val[2] =  -3.0;
  val[3] =   1.5;
  val[4] =   1.0;
  val[5] =   2.0;
  val[6] =   3.0;
  val[7] = 100.0;
  Array<double> valCopy = val; // deep copy

  const int expNumEntries = 5;
  const int    indExp[] = { 0,    1,   3,  -1,     0  };
  const double valExp[] = {42.0, -7.0, 1.5, 6.0, 100.0};

  // Test merge2 with default merge policy (add).
  {
    Array<int>::iterator indEnd = ind.end ();
    Array<double>::iterator valEnd = val.end ();
    merge2 (indEnd, valEnd, ind.begin (), indEnd, val.begin (), valEnd);

    const size_type newIndLen = indEnd - ind.begin ();
    const size_type newValLen = valEnd - val.begin ();

    TEST_EQUALITY( newIndLen, expNumEntries );
    TEST_EQUALITY( newValLen, expNumEntries );

    const bool indEq = std::equal (ind.begin (), indEnd, indExp);
    const bool valEq = std::equal (val.begin (), valEnd, valExp);

    TEST_EQUALITY( indEq, true );
    TEST_EQUALITY( valEq, true );

    if (! valEq) {
      out << "Input value range: " << toString (valCopy ()) << endl;
      out << "Expected output: "
          << toString (ArrayView<const double> ((const double*) valExp, expNumEntries))
          << endl;
      out << "Actual output: " << toString (val.view (0, newValLen)) << endl;
      Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
    }
  }

  ind = indCopy; // deep copy; restore original values
  val = valCopy; // deep copy; restore original values

  // Test merge2 with custom merge policy (also add).
  {
    Array<int>::iterator indEnd = ind.end ();
    Array<double>::iterator valEnd = val.end ();
    merge2 (indEnd, valEnd, ind.begin (), indEnd, val.begin (), valEnd, std::plus<double> ());

    const size_type newIndLen = indEnd - ind.begin ();
    const size_type newValLen = valEnd - val.begin ();

    TEST_EQUALITY( newIndLen, expNumEntries );
    TEST_EQUALITY( newValLen, expNumEntries );

    const bool indEq = std::equal (ind.begin (), indEnd, indExp);
    const bool valEq = std::equal (val.begin (), valEnd, valExp);

    TEST_EQUALITY( indEq, true );
    TEST_EQUALITY( valEq, true );

    if (! valEq) {
      out << "Input value range: " << toString (valCopy ()) << endl;
      out << "Expected output: "
          << toString (ArrayView<const double> ((const double*) valExp, expNumEntries))
          << endl;
      out << "Actual output: " << toString (val.view (0, newValLen)) << endl;
      Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
    }
  }
}


TEUCHOS_UNIT_TEST( MatrixMarket, CrsMatrixFileTest )
{
  success = testReadAndWriteFile (out, mapOutputFilename, matrixOutputFilename,
                                  mapInputFilename, matrixInputFilename);
}

