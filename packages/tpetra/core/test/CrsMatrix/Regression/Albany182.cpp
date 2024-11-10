// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include <algorithm>
#include <type_traits>
#include <vector>

// Test with 4 MPI processes.
// Test Export from the following overlapping CrsMatrix:
//
// PID  Row GID  Column GIDs
//   0        0  [0, 1,    3, 4, 5]
//   1        0      1, 2, 3, 5]
//   2        1  [3, 4,    6, 7, 8]
//   3        1  [   4, 5, 6,    8, 9]
//
// to the following nonoverlapping CrsMatrix:
//
// PID  Row GID  Column GIDs
//   0        0  [0, 1, 2, 3, 4, 5]
//   1           []
//   2        1  [3, 4, 5, 6, 7, 8, 9]
//   3           []
//
// Exercise the following target matrix cases:
//
//   - DynamicProfile, locally indexed
//   - Constant CrsGraph (locally indexed)
//
// The globally indexed use case is trickier, because the target
// matrix would construct its own column Map in that case.  The
// self-constructed target Map may differ from the column Map we use.
// That's OK, because that actually matches some application use cases
// (e.g., in Nalu, which constructs its own column Map that may not be
// the same as the Map that CrsGraph::makeColMap makes).  However, it
// complicates the test, because we have to convert the column indices
// to determine correctness.

namespace { // (anonymous)

  template<class ArrayType>
  void
  printArray (std::ostream& out, const ArrayType& x)
  {
    out << "[";
    const ptrdiff_t x_len = x.size ();
    for (ptrdiff_t k = 0; k < x_len; ++k) {
      out << x[k];
      if (k + ptrdiff_t (1) < x_len) {
        out << ",";
      }
    }
    out << "]";
  }

  template<class ArrayType1, class ArrayType2>
  void
  testArrayEquality (Teuchos::FancyOStream& out,
                     bool& success,
                     ArrayType1 x,
                     const size_t x_len,
                     const char x_name[],
                     ArrayType2 y,
                     const size_t y_len,
                     const char y_name[])
  {
    if (x_len != y_len) {
      TEST_ASSERT( false );
      out << "x_len = " << x_len << " != y_len = " << y_len << std::endl;
    }
    else {
      bool equal = true;
      for (size_t k = 0; k < x_len; ++k) {
        if (x[k] != y[k]) {
          equal = false;
          break;
        }
      }
      if (! equal) {
        TEST_ASSERT( false );
        out << "x = [";
        for (size_t k = 0; k < x_len; ++k) {
          out << x[k];
          if (k + size_t (1) < x_len) {
            out << ",";
          }
        }
        out << "] != y = [";
        for (size_t k = 0; k < y_len; ++k) {
          out << y[k];
          if (k + size_t (1) < y_len) {
            out << ",";
          }
        }
        out << "]" << std::endl;
      }
    }
  }

  template<class MapType>
  void
  printMapEntriesOnMyProcess (std::ostream& out,
                              const MapType& map)
  {
    using std::endl;
    typedef typename MapType::local_ordinal_type LO;
    typedef typename MapType::global_ordinal_type GO;

    if (map.getComm ().is_null ()) {
      return; // this process not participating in the Map
    }
    const int myRank = map.getComm ()->getRank ();
    out << "Proc " << myRank << ": [";
    const LO lclNumInds = static_cast<LO> (map.getLocalNumElements ());
    if (lclNumInds != 0) {
      for (LO lclInd = 0; lclInd < lclNumInds; ++lclInd) {
        const GO gblInd = map.getGlobalElement (lclInd);
        out << gblInd;
        if (lclInd + LO (1) < lclNumInds) {
          out << ",";
        }
      }
    }
    out << "]" << endl;
  }

  template<class MapType>
  void
  printMap (std::ostream& out,
            const MapType& map)
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = map.getComm ();
    if (! comm.is_null ()) {
      std::ostringstream os;
      printMapEntriesOnMyProcess (os, map);
      Tpetra::Details::gathervPrint (out, os.str (), *comm);
    }
  }

  template<class CrsMatrixType>
  void
  printMatrixEntriesOnMyProcess (std::ostream& out,
                                 const CrsMatrixType& A)
  {
    using std::endl;
    typedef typename CrsMatrixType::map_type map_type;
    typedef typename CrsMatrixType::local_ordinal_type LO;
    typedef typename CrsMatrixType::global_ordinal_type GO;

    const map_type& rowMap = * (A.getRowMap ());
    const int myRank = rowMap.getComm ()->getRank ();
    out << "Proc " << myRank << ": {";

    const LO lclNumRows = static_cast<LO> (rowMap.getLocalNumElements ());
    if (lclNumRows != 0) {
      if (A.isLocallyIndexed ()) {
        Teuchos::Array<GO> gblColInds;
        const map_type& colMap = * (A.getColMap ());
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const GO gblRow = rowMap.getGlobalElement (lclRow);
          out << "gblRow: " << gblRow;
          typename CrsMatrixType::local_inds_host_view_type lclColInds;
          typename CrsMatrixType::values_host_view_type vals;
          A.getLocalRowView (lclRow, lclColInds, vals);
          out << ": {lclCols: ";
          printArray (out, lclColInds);
          if (size_t(gblColInds.size ()) < lclColInds.size ()) {
            gblColInds.resize (lclColInds.size ());
          }
          for (size_t k = 0; k < lclColInds.size (); ++k) {
            gblColInds[k] = colMap.getGlobalElement (lclColInds[k]);
          }
          out << ", gblCols: ";
          printArray (out, gblColInds (0, lclColInds.size ()));
          out << ", vals: ";
          printArray (out, vals);
          out << "}";
          if (lclRow + LO (1) < lclNumRows) {
            out << ", ";
          }
        }
      }
      else {
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const GO gblRow = rowMap.getGlobalElement (lclRow);
          out << "gblRow: " << gblRow;
          typename CrsMatrixType::global_inds_host_view_type gblColInds;
          typename CrsMatrixType::values_host_view_type vals;
          A.getGlobalRowView (gblRow, gblColInds, vals);
          out << ": {gblCols: ";
          printArray (out, gblColInds);
          out << ", vals: ";
          printArray (out, vals);
          out << "}";
          if (lclRow + LO (1) < lclNumRows) {
            out << ", ";
          }
        }
      }
    }
    out << "}" << endl;
  }

  template<class MapType>
  Teuchos::RCP<const MapType>
  buildDomainMap (Teuchos::FancyOStream& out,
                  bool& success,
                  const bool verbose,
                  const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    typedef Tpetra::global_size_t GST;
    typedef MapType map_type;
    //typedef typename map_type::local_ordinal_type LO;
    typedef typename map_type::global_ordinal_type GO;

    if (comm.is_null ()) {
      out << "buildDomainMap: Input communicator is null!" << std::endl;
      TEST_ASSERT( false );
      return Teuchos::null;
    }
    else {
      const int myRank = comm->getRank ();
      const GO indexBase = 0;
      const GST gblNumInds_domain = 10;
      // [[0,1,2], [3,4,5], [6,7,8], [9]]
      const size_t lclNumInds_domain = (myRank >= 0 && myRank < 3) ?
        size_t (3) : (myRank == 3 ? size_t (1) : size_t (0));
      RCP<const map_type> domainMap (new map_type (gblNumInds_domain,
                                                   lclNumInds_domain,
                                                   indexBase, comm));
      if (verbose) {
        if (myRank == 0) {
          out << "Domain Map:" << std::endl;
        }
        Teuchos::OSTab tab1 (out);
        printMap (out, *domainMap);
      }
      return domainMap;
    }
  }

  template<class MapType>
  Teuchos::RCP<const MapType>
  buildRangeMap (Teuchos::FancyOStream& out,
                 bool& success,
                 const bool verbose,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    typedef Tpetra::global_size_t GST;
    typedef MapType map_type;
    //typedef typename map_type::local_ordinal_type LO;
    typedef typename map_type::global_ordinal_type GO;

    if (comm.is_null ()) {
      out << "buildRangeMap: Input communicator is null!" << std::endl;
      TEST_ASSERT( false );
      return Teuchos::null;
    }
    else {
      const int myRank = comm->getRank ();
      const GO indexBase = 0;
      const GST gblNumInds_range = 2;
      // [[0], [], [1], []]
      const size_t lclNumInds_range =
        (myRank == 0 || myRank == 2) ? size_t (1) : size_t (0);
      RCP<const map_type> rangeMap (new map_type (gblNumInds_range,
                                                  lclNumInds_range,
                                                  indexBase, comm));
      if (verbose) {
        if (myRank == 0) {
          out << "Range Map:" << std::endl;
        }
        Teuchos::OSTab tab1 (out);
        printMap (out, *rangeMap);
      }
      return rangeMap;
    }
  }

  template<class MapType>
  Teuchos::RCP<const MapType>
  buildOverlappingRowMap (Teuchos::FancyOStream& out,
                          bool& success,
                          const bool verbose,
                          const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    typedef Tpetra::global_size_t GST;
    typedef MapType map_type;
    typedef typename map_type::local_ordinal_type LO;
    typedef typename map_type::global_ordinal_type GO;

    if (comm.is_null ()) {
      out << "buildOverlappingRowMap: Input communicator is null!"
          << std::endl;
      TEST_ASSERT( false );
      return Teuchos::null;
    }
    else {
      const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
      std::vector<GO> rowMapInds;
      const int myRank = comm->getRank ();
      if (myRank == 0 || myRank == 1) {
        rowMapInds = std::vector<GO> {0};
      }
      else if (myRank == 2 || myRank == 3) {
        rowMapInds = std::vector<GO> {1};
      }
      const LO lclNumRowMapInds = static_cast<LO> (rowMapInds.size ());
      const GO indexBase = 0;
      RCP<const map_type> rowMap (new map_type (INVALID, rowMapInds.data (),
                                                lclNumRowMapInds,
                                                indexBase, comm));
      if (verbose) {
        if (myRank == 0) {
          out << "Overlapping row Map:" << std::endl;
        }
        Teuchos::OSTab tab1 (out);
        printMap (out, *rowMap);
      }
      return rowMap;
    }
  }

  template<class CrsMatrixType>
  void
  insertIntoOverlappingCrsMatrix (CrsMatrixType& A)
  {
    typedef typename CrsMatrixType::scalar_type ST;
    typedef typename CrsMatrixType::local_ordinal_type LO;
    typedef typename CrsMatrixType::global_ordinal_type GO;
    //typedef typename CrsMatrixType::map_type map_type;
    static_assert (std::is_same<ST, double>::value,
                   "CrsMatrixType::scalar_type must be double "
                   "in order for this test to work.");

    const int myRank = A.getMap ()->getComm ()->getRank ();
    if (myRank == 0) {
      const GO gblRow = 0;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {0, 1, 3, 4, 5};
      const ST vals[numEnt] = {0.0, 1.0, 3.0, 4.0, 5.0};
      A.insertGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 1) {
      const GO gblRow = 0;
      constexpr LO numEnt = 4;
      const GO inds[numEnt] = {1, 2, 3, 5};
      const ST vals[numEnt] = {1.0, 2.0, 3.0, 5.0};
      A.insertGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 2) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {3, 4, 6, 7, 8};
      const ST vals[numEnt] = {3.0, 4.0, 6.0, 7.0, 8.0};
      A.insertGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 3) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {4, 5, 6, 8, 9};
      const ST vals[numEnt] = {4.0, 5.0, 6.0, 8.0, 9.0};
      A.insertGlobalValues (gblRow, numEnt, vals, inds);
    }
  }

  template<class CrsGraphType>
  void
  insertIntoOverlappingCrsGraph (CrsGraphType& G)
  {
    typedef typename CrsGraphType::local_ordinal_type LO;
    typedef typename CrsGraphType::global_ordinal_type GO;

    const int myRank = G.getMap ()->getComm ()->getRank ();
    if (myRank == 0) {
      const GO gblRow = 0;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {0, 1, 3, 4, 5};
      G.insertGlobalIndices (gblRow, numEnt, inds);
    }
    else if (myRank == 1) {
      const GO gblRow = 0;
      constexpr LO numEnt = 4;
      const GO inds[numEnt] = {1, 2, 3, 5};
      G.insertGlobalIndices (gblRow, numEnt, inds);
    }
    else if (myRank == 2) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {3, 4, 6, 7, 8};
      G.insertGlobalIndices (gblRow, numEnt, inds);
    }
    else if (myRank == 3) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {4, 5, 6, 8, 9};
      G.insertGlobalIndices (gblRow, numEnt, inds);
    }
  }

  template<class CrsMatrixType>
  void
  fillIntoOverlappingCrsMatrix (CrsMatrixType& A)
  {
    typedef typename CrsMatrixType::scalar_type ST;
    typedef typename CrsMatrixType::local_ordinal_type LO;
    typedef typename CrsMatrixType::global_ordinal_type GO;
    //typedef typename CrsMatrixType::map_type map_type;
    static_assert (std::is_same<ST, double>::value,
                   "CrsMatrixType::scalar_type must be double "
                   "in order for this test to work.");

    const int myRank = A.getMap ()->getComm ()->getRank ();
    if (myRank == 0) {
      const GO gblRow = 0;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {0, 1, 3, 4, 5};
      const ST vals[numEnt] = {0.0, 1.0, 3.0, 4.0, 5.0};
      (void) A.replaceGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 1) {
      const GO gblRow = 0;
      constexpr LO numEnt = 4;
      const GO inds[numEnt] = {1, 2, 3, 5};
      const ST vals[numEnt] = {1.0, 2.0, 3.0, 5.0};
      (void) A.replaceGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 2) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {3, 4, 6, 7, 8};
      const ST vals[numEnt] = {3.0, 4.0, 6.0, 7.0, 8.0};
      (void) A.replaceGlobalValues (gblRow, numEnt, vals, inds);
    }
    else if (myRank == 3) {
      const GO gblRow = 1;
      constexpr LO numEnt = 5;
      const GO inds[numEnt] = {4, 5, 6, 8, 9};
      const ST vals[numEnt] = {4.0, 5.0, 6.0, 8.0, 9.0};
      (void) A.replaceGlobalValues (gblRow, numEnt, vals, inds);
    }
  }

  template<class CrsMatrixType>
  void
  testEntriesOfOverlappingCrsMatrix (Teuchos::FancyOStream& out,
                                     bool& success,
                                     const bool /* verbose */,
                                     const CrsMatrixType& A)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef typename CrsMatrixType::scalar_type ST;
    typedef typename CrsMatrixType::local_ordinal_type LO;
    typedef typename CrsMatrixType::global_ordinal_type GO;
    typedef typename CrsMatrixType::map_type map_type;
    static_assert (std::is_same<ST, double>::value,
                   "CrsMatrixType::scalar_type must be double "
                   "in order for this test to work.");

    out << "Test entries of overlapping CrsMatrix" << std::endl;
    Teuchos::OSTab tab1 (out);

    // mfh 28 Sep 2017: Historically, an Epetra_CrsMatrix or a
    // Tpetra::CrsMatrix could be neither globally nor locally indexed
    // on a process, if it had no entries on that process.
    const bool globallyIndexed = A.isGloballyIndexed ();
    const bool locallyIndexed = A.isLocallyIndexed ();
    const map_type& rowMap = * (A.getRowMap ());
    RCP<const map_type> colMap =
      locallyIndexed ? A.getColMap () : Teuchos::null;
    RCP<const Comm<int> > comm = A.getMap ()->getComm ();
    const int myRank = comm->getRank ();

    // These four things will get set below; they are process-specific
    GO gblRow = 0;
    LO numEnt_expected = 0;
    std::vector<GO> gblInds_expected;
    std::vector<ST> vals_expected;

    if (myRank == 0) {
      gblRow = 0;
      numEnt_expected = 5;
      gblInds_expected = std::vector<GO> {0, 1, 3, 4, 5};
      vals_expected = std::vector<ST> {0.0, 1.0, 3.0, 4.0, 5.0};
    }
    else if (myRank == 1) {
      gblRow = 0;
      numEnt_expected = 4;
      gblInds_expected = std::vector<GO> {1, 2, 3, 5};
      vals_expected = std::vector<ST> {1.0, 2.0, 3.0, 5.0};
    }
    else if (myRank == 2) {
      gblRow = 1;
      numEnt_expected = 5;
      gblInds_expected = std::vector<GO> {3, 4, 6, 7, 8};
      vals_expected = std::vector<ST> {3.0, 4.0, 6.0, 7.0, 8.0};
    }
    else if (myRank == 3) {
      gblRow = 1;
      numEnt_expected = 5;
      gblInds_expected = std::vector<GO> {4, 5, 6, 8, 9};
      vals_expected = std::vector<ST> {4.0, 5.0, 6.0, 8.0, 9.0};
    }

    if (myRank >= 0 && myRank < 4) {
      if (locallyIndexed) {
        out << "Overlapping CrsMatrix is locally indexed" << endl;
        Teuchos::OSTab tab2 (out);

        const GO lclRow = rowMap.getLocalElement (gblRow);
        typename CrsMatrixType::crs_graph_type::local_inds_host_view_type
                 lclInds;
        typename CrsMatrixType::values_host_view_type vals;
        A.getLocalRowView(lclRow, lclInds, vals);
        LO numEnt = lclInds.extent(0);

        TEST_EQUALITY( numEnt, numEnt_expected );
        if (numEnt == numEnt_expected) {
          // Convert to global column indices for comparison.  Note
          // that the overlapping matrix's column Map is not the same
          // as the nonoverlapping matrix's column Map, so we can't
          // compare local column indices.
          std::vector<GO> gblInds_got (numEnt);
          for (LO k = 0; k < numEnt; ++k) {
            gblInds_got[k] = colMap->getGlobalElement (lclInds[k]);
          }
          // Sort the global column indices, since we don't want to
          // rely on a promise of ordering.  Apply the same
          // permutation to the values, so we can compare.
          std::vector<ST> vals_got(numEnt);
          for (LO k = 0; k < numEnt; ++k) vals_got[k] = vals[k];
          Tpetra::sort2 (gblInds_got.begin (), gblInds_got.end (),
                         vals_got.begin ());
          testArrayEquality (out, success,
                             gblInds_got.data (), gblInds_got.size (),
                             "gblInds_got",
                             gblInds_expected.data (), gblInds_expected.size (),
                             "gblInds_expected");
          testArrayEquality (out, success,
                             vals_got.data (), vals_got.size (),
                             "vals_got",
                             vals_expected.data (), vals_expected.size (),
                             "vals_expected");
        }
      }
      else if (globallyIndexed) {
        out << "Overlapping CrsMatrix is globally indexed" << endl;
        Teuchos::OSTab tab2 (out);

        typename CrsMatrixType::global_inds_host_view_type gblInds;
        typename CrsMatrixType::values_host_view_type vals;
        A.getGlobalRowView (gblRow, gblInds, vals);
        const LO numEnt = static_cast<LO> (gblInds.size ());

        TEST_EQUALITY( numEnt, numEnt_expected );
        if (numEnt == numEnt_expected) {
          // Copy the input matrix data, and sort the global column
          // indices, since we don't want to rely on a promise of
          // ordering.  Apply the same permutation to the values, so
          // we can compare.
          std::vector<GO> gblInds_got(numEnt);
          std::vector<ST> vals_got(numEnt);
          for (LO k = 0; k < numEnt; ++k) {
            gblInds_got[k] = gblInds[k];
            vals_got[k] = vals[k];
          }
          Tpetra::sort2 (gblInds_got.begin (), gblInds_got.end (),
                         vals_got.begin ());
          testArrayEquality (out, success, gblInds_got.data (),
                             gblInds_got.size (), "gblInds_got",
                             gblInds_expected.data (),
                             gblInds_expected.size (),
                             "gblInds_expected");
          testArrayEquality (out, success, vals_got.data (),
                             vals_got.size (), "vals_got",
                             vals_expected.data (),
                             vals_expected.size (),
                             "vals_expected");
        }
      }
      else {
        out << "Overlapping CrsMatrix is neither locally nor globally "
          "indexed; that means the matrix has no entries on this process, "
          "which should not happen!" << endl;
        Teuchos::OSTab tab2 (out);
        TEST_ASSERT( false );
      }
    }

    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (myRank == 0) {
        out << "FAILED to build overlapping CrsMatrix correctly!" << std::endl;
      }
    }
  }

  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  buildOverlappingCrsMatrix (Teuchos::FancyOStream& out,
                             bool& success,
                             const bool staticGraph,
                             const bool verbose,
                             const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    //typedef Tpetra::global_size_t GST;
    //typedef typename CrsMatrixType::local_ordinal_type LO;
    //typedef typename CrsMatrixType::global_ordinal_type GO;
    typedef typename CrsMatrixType::map_type map_type;

    RCP<const map_type> rangeMap = buildRangeMap<map_type> (out, success, verbose, comm);
    RCP<const map_type> domainMap = buildDomainMap<map_type> (out, success, verbose, comm);
    RCP<const map_type> rowMap = buildOverlappingRowMap<map_type> (out, success, verbose, comm);

    RCP<CrsMatrixType> A;
    {
      using Teuchos::rcp_const_cast;
      typedef typename CrsMatrixType::crs_graph_type crs_graph_type;
      RCP<crs_graph_type> G (new crs_graph_type (rowMap, 5));
      insertIntoOverlappingCrsGraph (*G);
      G->fillComplete (domainMap, rangeMap);
      // typedef typename CrsMatrixType::device_type::execution_space execution_space;
      // execution_space().fence ();
      A = rcp (new CrsMatrixType (rcp_const_cast<const crs_graph_type> (G)));
      fillIntoOverlappingCrsMatrix (*A);
    }
    A->fillComplete (domainMap, rangeMap);
    // typedef typename CrsMatrixType::device_type::execution_space execution_space;
    // execution_space().fence ();
    testEntriesOfOverlappingCrsMatrix (out, success, verbose, *A);

    return A;
  }

  template<class MapType>
  Teuchos::RCP<const MapType>
  expectedColumnMapFromNonoverlappingCrsMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    typedef Tpetra::global_size_t GST;
    typedef MapType map_type;
    typedef typename map_type::local_ordinal_type LO;
    typedef typename map_type::global_ordinal_type GO;

    // Expected column Map:
    // [[0,1,2,3,4,5], [], [3,4,5,6,7,8,9], []].

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    std::vector<GO> inds;
    const int myRank = comm->getRank ();
    if (myRank == 0) {
      inds = std::vector<GO> {0, 1, 2, 3, 4, 5};
    }
    else if (myRank == 2) {
      inds = std::vector<GO> {3, 4, 5, 6, 7, 8, 9};
    }
    const LO lclNumInds = static_cast<LO> (inds.size ());
    const GO indexBase = 0;
    RCP<const map_type> map (new map_type (INVALID, inds.data (), lclNumInds,
                                           indexBase, comm));
    return map;
  }

  template<class GO>
  std::pair<std::vector<double>, std::vector<GO> >
  expectedGblColIndsAndValsInNonoverlappingCrsMatrix (const int myRank)
  {
    if (myRank == 0) {
      return std::make_pair (std::vector<double> {0.0, 2.0, 2.0, 6.0, 4.0, 10.0},
                             std::vector<GO>     {0,   1,   2,   3,   4,    5});
    }
    else if (myRank == 2) {
      return std::make_pair (std::vector<double> {3.0, 8.0, 5.0, 12.0, 7.0, 16.0, 9.0},
                             std::vector<GO>     {3,   4,   5,    6,   7,    8,   9});
    }
    else {
      return std::make_pair (std::vector<double> {}, std::vector<GO> {});
    }
  }

  template<class MapType>
  bool
  testNonoverlappingColumnMapIsAsExpected (Teuchos::FancyOStream& out,
                                           bool& success,
                                           const MapType* const colMapPtr,
                                           const MapType& colMap_expected)
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;

    const Comm<int>& comm = * (colMap_expected.getComm ());
    const int myRank = comm.getRank ();
    const int colMapNonNull_lcl = (colMapPtr == NULL) ? 0 : 1;
    int colMapNonNull_gbl = 0; // output argument
    reduceAll<int, int> (comm, REDUCE_MIN, colMapNonNull_lcl,
                         outArg (colMapNonNull_gbl));
    TEST_EQUALITY_CONST( colMapNonNull_gbl, 1 );
    if (colMapNonNull_gbl == 1) {
      const bool colMapsSame = colMap_expected.isSameAs (*colMapPtr);
      TEST_ASSERT( colMapsSame );
      if (! colMapsSame) {
        if (myRank == 0) {
          out << "The expected column Map and the nonoverlapping "
            "matrix's column Map are not the same!" << endl;
          out << "Expected column Map:" << endl;
        }
        printMap (out, colMap_expected);
        if (myRank == 0) {
          out << "Actual column Map:" << endl;
        }
        printMap (out, *colMapPtr);
      }
      return colMapsSame;
    }
    else {
      out << "The matrix's column Map is null on one or more processes!"
          << endl;
      return false;
    }
  }

  template<class CrsMatrixType>
  void
  testEntriesOfNonoverlappingCrsMatrix (Teuchos::FancyOStream& out,
                                        bool& success,
                                        const CrsMatrixType& A_nonoverlapping,
                                        const typename CrsMatrixType::map_type& colMap_expected)
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    //typedef typename CrsMatrixType::map_type map_type;
    typedef typename CrsMatrixType::local_ordinal_type LO;
    typedef typename CrsMatrixType::global_ordinal_type GO;

    out << "Test entries of nonoverlapping CrsMatrix" << endl;
    Teuchos::OSTab tab1 (out);

    const Comm<int>& comm = * (A_nonoverlapping.getComm ());
    const int myRank = comm.getRank ();
    // We know at this point that the matrix is locally indexed,
    // since it has been fill-completed (at least) once.
    //const map_type& colMap = * (A_nonoverlapping.getColMap ());

    std::vector<double> vals_expected;
    std::vector<LO> lclColInds_expected;
    std::vector<GO> gblColInds_expected;
    LO numEnt_expected = 0;
    {
      std::pair<std::vector<double>, std::vector<GO> > expectedData =
        expectedGblColIndsAndValsInNonoverlappingCrsMatrix<GO> (myRank);
      std::swap (vals_expected, expectedData.first);
      std::swap (gblColInds_expected, expectedData.second);
      TEST_EQUALITY( vals_expected.size (), gblColInds_expected.size () );
      numEnt_expected = static_cast<LO> (vals_expected.size ());

      // Convert the expected global column indices to local column
      // indices.  Use the expected column Map, not the actual one,
      // just in case the two column Maps differ.  (That would
      // indicate a bug, but we still want to extract useful
      // information in that case.)
      lclColInds_expected = std::vector<LO> (numEnt_expected);
      for (LO k = 0; k < numEnt_expected; ++k) {
        const LO gblColInd = gblColInds_expected[k];
        lclColInds_expected[k] = colMap_expected.getLocalElement (gblColInd);
      }

      // Sort the expected local column indices, and apply the sort
      // permutation to the expected values and the corresponding
      // global column indices.  We'll need that because the matrix
      // rows are sorted by local column index, not by global column
      // index.
      Tpetra::sort3 (lclColInds_expected.begin (),
                     lclColInds_expected.end (),
                     vals_expected.begin (),
                     gblColInds_expected.begin ());
    }

    if (myRank == 0 || myRank == 2) {
      const GO gblRow = (myRank == 0) ? GO (0) : GO (1);
      const LO lclRow =
        A_nonoverlapping.getRowMap ()->getLocalElement (gblRow);
      typename CrsMatrixType::local_inds_host_view_type lclColInds;
      typename CrsMatrixType::values_host_view_type vals;
      A_nonoverlapping.getLocalRowView(lclRow, lclColInds, vals);
      LO numEnt = lclColInds.extent(0);
      // Sort the input matrix's row data, so that this test does
      // not depend on any assumption of sorted rows.
      std::vector<LO> lclColInds_got (numEnt);
      std::vector<double> vals_got(numEnt);
      for (LO k = 0; k < numEnt; ++k) {
        lclColInds_got[k] = lclColInds[k];
        vals_got[k] = vals[k];
      }
      Tpetra::sort2 (lclColInds_got.begin (), lclColInds_got.end (),
                     vals_got.begin ());
      testArrayEquality (out, success, lclColInds_got.data (),
                         lclColInds_got.size (), "lclColInds_got",
                         lclColInds_expected.data (),
                         lclColInds_expected.size (),
                         "lclColInds_expected");
      testArrayEquality (out, success, vals_got.data (), vals_got.size (),
                         "vals_got", vals_expected.data (),
                         vals_expected.size (), "vals_expected");
    }
    else {
      const LO lclNumRows =
        static_cast<LO> (A_nonoverlapping.getLocalNumRows ());
      TEST_EQUALITY_CONST( lclNumRows, LO (0) );
    }

    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (gblSuccess != 1) {
      if (myRank == 0) {
        out << "Overlapping matrix is wrong!  Here is what it should be:"
            << endl;
      }
      {
        std::ostringstream os;
        os << "Proc " << myRank << ": {";
        if (myRank == 0 || myRank == 2) {
          const GO gblRow = (myRank == 0) ? GO (0) : GO (1);
          os << "gblRow: " << gblRow << ": {lclCols: ";
          printArray (os, lclColInds_expected);
          os << ", gblCols: ";
          printArray (os, gblColInds_expected);
          os << ", vals: ";
          printArray (os, vals_expected);
        }
        os << "}" << endl;
        Tpetra::Details::gathervPrint (out, os.str (), comm);
      }
      if (myRank == 0) {
        out << "Here is the actual overlapping matrix:" << endl;
      }
      {
        std::ostringstream os;
        printMatrixEntriesOnMyProcess (os, A_nonoverlapping);
        Tpetra::Details::gathervPrint (out, os.str (), comm);
      }
    }
  }

  template<class CrsMatrixType>
  void
  createAndExportAndTestCrsMatrix (Teuchos::FancyOStream& out,
                                   bool& success,
                                   const CrsMatrixType& A_overlapping)
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef typename CrsMatrixType::map_type map_type;
    typedef typename map_type::local_ordinal_type LO;
    typedef typename map_type::global_ordinal_type GO;
    //typedef typename map_type::device_type::execution_space execution_space;
    typedef Tpetra::Export<typename map_type::local_ordinal_type,
      typename map_type::global_ordinal_type,
      typename map_type::node_type> export_type;
    int lclSuccess = 1; // to be modified below
    int gblSuccess = 0; // output argument

    RCP<const Comm<int> > comm = A_overlapping.getMap ()->getComm ();
    const int myRank = comm->getRank ();

    Teuchos::OSTab tab0 (out);
    out << "createAndExportAndTestCrsMatrix" << endl;
    Teuchos::OSTab tab1 (out);

    TEST_ASSERT( A_overlapping.isFillComplete () );
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (myRank == 0) {
        out << "Test FAILED on one or more processes; exiting early!" << endl;
      }
      return;
    }

    RCP<const map_type> domMap = A_overlapping.getDomainMap ();
    RCP<const map_type> ranMap = A_overlapping.getRangeMap ();
    RCP<const map_type> rowMap_nonoverlapping = ranMap;

    RCP<const map_type> colMap_expected =
      expectedColumnMapFromNonoverlappingCrsMatrix<map_type> (comm);

    if (false) {
      out << "Target matrix is {DynamicProfile, globally indexed}" << endl;
      Teuchos::OSTab tab2 (out);

      out << "Create target matrix (A_nonoverlapping)" << endl;
      RCP<CrsMatrixType> A_nonoverlapping =
        rcp (new CrsMatrixType (rowMap_nonoverlapping, 0));

      out << "Create Export object" << endl;
      export_type exp (A_overlapping.getRowMap (), rowMap_nonoverlapping);

      out << "Export from source matrix to target matrix" << endl;
      A_nonoverlapping->doExport (A_overlapping, exp, Tpetra::ADD);

      out << "Call fillComplete on the target matrix" << endl;
      A_nonoverlapping->fillComplete (domMap, ranMap);
      //execution_space().fence ();

      out << "Test target matrix's column Map" << endl;
      const map_type* colMapPtr = A_nonoverlapping->getColMap ().getRawPtr ();
      const bool colMapsSame =
        testNonoverlappingColumnMapIsAsExpected (out, success, colMapPtr,
                                                 *colMap_expected);
      if (! colMapsSame && myRank == 0) {
        out << "Column Maps are not as expected, so it's doubtful that the "
          "nonoverlapping matrix is right.  Nevertheless, we'll continue, "
          "because looking at the actual matrix entries could help with "
          "debugging." << endl;
      }
      out << "Test target matrix's entries" << endl;
      testEntriesOfNonoverlappingCrsMatrix (out, success,
                                            *A_nonoverlapping,
                                            *colMap_expected);
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        if (myRank == 0) {
          out << "Test FAILED on one or more processes; exiting early!" << endl;
        }
        return;
      }
      out << "Target matrix is correct!" << endl;
    }

    {
      out << ">>> Target matrix is locally indexed" << endl;
      Teuchos::OSTab tab2 (out);

      const size_t maxNumEntPerRow = 10; // needs to be an upper bound
      RCP<CrsMatrixType> A_nonoverlapping =
        rcp (new CrsMatrixType (rowMap_nonoverlapping, colMap_expected,
                                maxNumEntPerRow));
      export_type exp (A_overlapping.getRowMap (), rowMap_nonoverlapping);
      A_nonoverlapping->doExport (A_overlapping, exp, Tpetra::ADD);
      A_nonoverlapping->fillComplete (domMap, ranMap);
      //execution_space().fence ();
      const map_type* colMapPtr = A_nonoverlapping->getColMap ().getRawPtr ();
      const bool colMapsSame =
        testNonoverlappingColumnMapIsAsExpected (out, success, colMapPtr,
                                                 *colMap_expected);
      if (! colMapsSame && myRank == 0) {
        out << "Column Maps are not as expected, so it's doubtful that the "
          "nonoverlapping matrix is right.  Nevertheless, we'll continue, "
          "because looking at the actual matrix entries could help with "
          "debugging." << endl;
      }
      testEntriesOfNonoverlappingCrsMatrix (out, success,
                                            *A_nonoverlapping,
                                            *colMap_expected);
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        out << "Test FAILED on one or more processes; exiting early!" << endl;
        return;
      }
    }

    {
      out << ">>> Target matrix has a const CrsGraph, "
        "and has been fill-completed once" << endl;
      Teuchos::OSTab tab2 (out);

      typedef typename CrsMatrixType::crs_graph_type crs_graph_type;
      const size_t maxNumEntPerRow = (myRank == 0) ? size_t (6) :
        ((myRank == 2) ? size_t (7) : size_t (0));
      RCP<crs_graph_type> G_nonoverlapping =
        rcp (new crs_graph_type (rowMap_nonoverlapping, colMap_expected,
                                 maxNumEntPerRow));
      if (myRank == 0) {
        const GO gblRow = 0;
        std::vector<GO> gblColInds {0, 1, 2, 3, 4, 5};
        const LO numEnt = static_cast<LO> (gblColInds.size ());
        G_nonoverlapping->insertGlobalIndices (gblRow, numEnt,
                                               gblColInds.data ());
      }
      else if (myRank == 2) {
        const GO gblRow = 1;
        std::vector<GO> gblColInds {3, 4, 5, 6, 7, 8, 9};
        const LO numEnt = static_cast<LO> (gblColInds.size ());
        G_nonoverlapping->insertGlobalIndices (gblRow, numEnt,
                                               gblColInds.data ());
      }
      //execution_space().fence ();
      G_nonoverlapping->fillComplete (domMap, ranMap);
      const map_type* colMapPtr_G = G_nonoverlapping->getColMap ().getRawPtr ();
      const bool gotGraphColMapRight =
        testNonoverlappingColumnMapIsAsExpected (out, success, colMapPtr_G,
                                                 *colMap_expected);
      if (! gotGraphColMapRight && myRank == 0) {
        out << "Hm, the CrsGraph G_nonoverlapping has the wrong column Map.  "
          "This almost certainly means that the test is totally wrong.  "
          "We'll continue, in hopes of learning something." << endl;
      }
      RCP<CrsMatrixType> A_nonoverlapping (new CrsMatrixType (G_nonoverlapping));
      export_type exp (A_overlapping.getRowMap (), rowMap_nonoverlapping);
      A_nonoverlapping->doExport (A_overlapping, exp, Tpetra::ADD);
      A_nonoverlapping->fillComplete (domMap, ranMap);
      A_nonoverlapping->resumeFill ();
      //execution_space().fence ();
      const map_type* colMapPtr = A_nonoverlapping->getColMap ().getRawPtr ();
      const bool colMapsSame =
        testNonoverlappingColumnMapIsAsExpected (out, success, colMapPtr,
                                                 *colMap_expected);
      if (! colMapsSame && myRank == 0) {
        out << "Column Maps are not as expected, so it's doubtful that the "
          "nonoverlapping matrix is right.  Nevertheless, we'll continue, "
          "because looking at the actual matrix entries could help with "
          "debugging." << endl;
      }
      testEntriesOfNonoverlappingCrsMatrix (out, success,
                                            *A_nonoverlapping,
                                            *colMap_expected);
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        out << "Test FAILED on one or more processes; exiting early!" << endl;
        return;
      }
    }
  }

  template<class CrsMatrixType>
  void
  testCrsMatrixExport (Teuchos::FancyOStream& out,
                       bool& success,
                       const bool verbose,
                       const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    using Teuchos::RCP;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef CrsMatrixType crs_matrix_type;
    int lclSuccess = 1; // to be revised below
    int gblSuccess = 0; // output argument

    const int myRank = comm->getRank ();
    Teuchos::OSTab tab0 (out);
    out << "testCrsMatrixExport" << endl;
    Teuchos::OSTab tab1 (out);

    for (bool staticGraph : {true})
    {
      out << "Source matrix: staticGraph=" << (staticGraph ? "true" : "false")
          << endl;
      RCP<crs_matrix_type> A_overlapping =
        buildOverlappingCrsMatrix<crs_matrix_type> (out, success,
                                                    staticGraph,
                                                    verbose,
                                                    comm);
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        out << "FAILED to build overlapping CrsMatrix correctly; "
          "returning from test early!" << std::endl;
        return;
      }

      TEST_ASSERT( A_overlapping->isFillComplete () );
      TEST_EQUALITY( staticGraph, A_overlapping->isStaticGraph () );
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        out << "FAILED to set up overlapping CrsMatrix correctly; "
          "returning from test early!" << std::endl;
        return;
      }

      createAndExportAndTestCrsMatrix (out, success, *A_overlapping);
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        if (myRank == 0) {
          out << "A_overlapping (source matrix of the Export) staticGraph="
              << (staticGraph ? "true" : "false") << " case FAILED!" << endl;
        }
      }
    }

  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, Albany182, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsMatrix<double, LO, GO, Node> crs_matrix_type;

    RCP<const Comm<int> > comm = Tpetra::TestingUtilities::getDefaultComm ();
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    testCrsMatrixExport<crs_matrix_type> (out, success, verbose, comm);
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix, Albany182, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  //TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

  typedef ::Tpetra::Map<>::node_type default_node_type;
#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
  UNIT_TEST_GROUP( int, longlong, default_node_type )
#else
  typedef ::Tpetra::Map<>::local_ordinal_type default_local_ordinal_type;
  typedef ::Tpetra::Map<>::global_ordinal_type default_global_ordinal_type;
  UNIT_TEST_GROUP( default_local_ordinal_type, default_global_ordinal_type, default_node_type )
#endif // HAVE_TPETRA_INST_INT_LONG_LONG

} // namespace (anonymous)
