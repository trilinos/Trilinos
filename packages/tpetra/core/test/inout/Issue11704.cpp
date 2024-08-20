// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_ETIHelperMacros.h"
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Util.hpp> // sort2, merge2

/*! \file Issue11704.cpp

   https://github.com/trilinos/Trilinos/issues/11704

    Issue 11704 says that the crs reader only populates the first and last
   entries of the rowPtr array, and then also any entry up to the last non-zero.
    This means that if final rows of a matrix are empty, the last row
    of the matrix will almost certainly have the wrong length (it will go from 0
   to nnz)

    Test this by creating some matrices with all non-zeros in the first row, and
   then asserting that it is true for the resulting matrix
*/

namespace { // anonymous

using std::endl;
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
using Tpetra::global_size_t;
using Tpetra::TestingUtilities::arcp_from_view;

struct Test {
  const char *contents;
  size_t numRows;
  size_t numCols;
  size_t nnz;
  size_t populatedRow;
};

const Test EMPTY = {"%%MatrixMarket matrix coordinate real general\n"
                    "5 5 0\n",
                    5, 5, 0, 0 /*don't care*/};

const Test FIRST = {"%%MatrixMarket matrix coordinate real general\n"
                    "5 5 13\n"
                    "1 1  2.0\n",
                    5, 5, 1, 0};

template <class ScalarType, class LocalOrdinalType, class GlobalOrdinalType,
          class NodeType>
bool testReadSparseStream(Teuchos::FancyOStream &out, const Test &test,
                          RCP<const Comm<int>> &comm) {

  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  ;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;
  using inds_type = typename crs_matrix_type::local_inds_host_view_type;
  using vals_type = typename crs_matrix_type::values_host_view_type;

  const bool callFillComplete = true;
  const bool tolerant = false;
  const bool debug = false;

  out << "Original sparse matrix:" << endl;
  out << test.contents << endl;
  std::istringstream iss(test.contents);

  out << "Creating the row Map" << endl;
  RCP<const map_type> rowMap =
      rcp(new map_type(test.numRows, 0, comm, Tpetra::GloballyDistributed));

  out << "Reading in the matrix" << endl;
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  RCP<crs_matrix_type> A =
      reader_type::readSparse(iss, rowMap, colMap, domainMap, rangeMap,
                              callFillComplete, tolerant, debug);

  bool success = true;
  TEUCHOS_TEST_EQUALITY(A->getGlobalNumEntries(), test.nnz, out, success);
  TEUCHOS_TEST_EQUALITY(A->getGlobalNumRows(), test.numRows, out, success);
  TEUCHOS_TEST_EQUALITY(A->getGlobalNumCols(), test.numCols, out, success);

  // use rank 0 to check that all non-zeros are in the correct rows
  for (size_t lrow = 0; lrow < A->getLocalNumRows(); ++lrow) {
    GlobalOrdinalType grow = rowMap->getGlobalElement(lrow);
    inds_type inds;
    vals_type vals;
    A->getLocalRowView(lrow, inds, vals);
    out << "test lrow=" << lrow << " grow=" << grow << std::endl;
    if (grow == GlobalOrdinalType(test.populatedRow)) {
      TEUCHOS_TEST_EQUALITY(inds.size(), test.nnz, out, success);
    } else {
      TEUCHOS_TEST_EQUALITY(inds.size(), size_t(0), out, success);
    }
  }

  return success;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool testReadSparse(Teuchos::FancyOStream &out) {

  out << "Test: https://github.com/trilinos/Trilinos/issues/11704" << std::endl;
  OSTab tab1(out);

  RCP<const Comm<int>> comm = Tpetra::getDefaultComm();

  bool success = true;
  success = success &&
            testReadSparseStream<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
                out, EMPTY, comm);
  success = success &&
            testReadSparseStream<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
                out, FIRST, comm);
  return success;
}

} // namespace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrixOutputInput, Issue11704, ST, LO, GO,
                                  NT) {
  success = testReadSparse<ST, LO, GO, NT>(out);
}

#if defined(HAVE_TPETRA_INST_DOUBLE)
#define UNIT_TEST_GROUP(LO, GO, NODE)                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixOutputInput, Issue11704,       \
                                       double, LO, GO, NODE)

#elif defined(HAVE_TPETRA_INST_FLOAT)
#define UNIT_TEST_GROUP(LO, GO, NODE)                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixOutputInput, Issue11704,       \
                                       float, LO, GO, NODE)                    \
#else
#define UNIT_TEST_GROUP(LO, GO, NODE)
#endif

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)
