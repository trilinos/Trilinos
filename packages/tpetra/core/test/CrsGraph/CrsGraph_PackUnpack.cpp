// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_packCrsGraph.hpp"
#include "Tpetra_Details_unpackCrsGraphAndCombine.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <random>
#include <set>

namespace { // anonymous

#define NUM_ROW_PER_PROC 4
#define NUM_NZ_COLS 100

using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::outArg;
using Tpetra::Details::gathervPrint;
using Tpetra::Details::packCrsGraph;
using std::endl;

template<class T>
bool
essentially_equal(T a, T b) {
  typedef Kokkos::ArithTraits<T> KAT;
  const auto eps = KAT::eps();
  return KAT::abs(a - b) <= ( (KAT::abs(a) > KAT::abs(b) ? KAT::abs(b) : KAT::abs(a)) * eps);
}

template<class CrsGraphType>
Teuchos::RCP<CrsGraphType>
generate_test_graph(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  typedef CrsGraphType crs_graph_type;
  typedef typename crs_graph_type::local_ordinal_type LO;
  typedef typename crs_graph_type::global_ordinal_type GO;
  typedef typename crs_graph_type::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> MapType;

  const int world_rank = comm->getRank();

  const LO num_row_per_proc = NUM_ROW_PER_PROC; // 4;
  Array<GO> row_gids(num_row_per_proc);

  GO start = static_cast<GO>(num_row_per_proc * world_rank);
  for (LO i=0; i<num_row_per_proc; ++i) {
    row_gids[i] = static_cast<GO>(start+i);
  }

  // Create random, unique column GIDs.
  const LO max_num_ent_per_row = NUM_NZ_COLS; //100;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_int_distribution<GO>  distr(1, 2000);
  std::set<GO> col_gids_set;
  typedef typename std::set<GO>::size_type SGO;
  SGO num_gids_in_set = static_cast<SGO>(max_num_ent_per_row);
  col_gids_set.insert(0);
  int num_passes = 0;
  while (col_gids_set.size() < num_gids_in_set && num_passes <= 2000) {
    col_gids_set.insert(distr(generator));
    num_passes += 1;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(col_gids_set.size() != num_gids_in_set,
      std::runtime_error, "Random column IDs not generated");

  Array<GO> col_gids(col_gids_set.begin(), col_gids_set.end());

  const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
  RCP<const MapType> row_map = rcp(new MapType(INVALID, row_gids(), 0, comm));
  RCP<const MapType> col_map = rcp(new MapType(INVALID, col_gids(), 0, comm));

  auto A = rcp(new crs_graph_type(row_map, col_map, max_num_ent_per_row));

  Array<LO> columns(max_num_ent_per_row);
  for (LO j=0; j<max_num_ent_per_row; ++j) columns[j] = j;
  for (LO i=0; i<num_row_per_proc; ++i) {
    //LO lcl_row = static_cast<LO>(start + i); // unused
    A->insertLocalIndices(i, columns());
  }
  A->fillComplete(col_map, row_map);

  return A;
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, PackThenUnpackAndCombine, LO, GO, NT)
{
  typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
  typedef typename crs_graph_type::packet_type packet_type;

  int lclSuccess = 1; // to be revised below
  int gblSuccess = 0; // output argument

  RCP<const Comm<int> > comm = getDefaultComm();
  const int world_rank = comm->getRank();

  out << "Creating graph" << endl;

  auto A = generate_test_graph<crs_graph_type> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Preparing arguments for packCrsGraph" << endl;

  LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
  Array<LO> exportLIDs (num_loc_rows); // input argument
  for (LO i=0; i < num_loc_rows; ++i) {
    exportLIDs[i] = static_cast<LO>(i); // pack all the rows
  }
  Array<packet_type> exports; // output argument; to be realloc'd
  Array<size_t> numPacketsPerLID (num_loc_rows, 0); // output argument
  size_t constantNumPackets; // output argument

  out << "Calling packCrsGraph" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsGraph<LO,GO,NT>(*A, exports, numPacketsPerLID(), exportLIDs(),
          constantNumPackets);
      local_op_ok = 1;
    } catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }
    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "packCrsGraph reported an error!" << endl;
      }
      gathervPrint (out, msg.str(), *comm);
      out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Now make sure that the pack is correct by creating an empty graph and
  // unpacking in to it.  The graph should end up being the same as the above graph.
  out << "Building second graph" << endl;
  RCP<crs_graph_type> B = rcp(new crs_graph_type(row_map, col_map, A->getLocalNumEntries()));

#if 0
  out << "Calling unpackCrsGraphAndCombine" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      unpackCrsGraphAndCombine<LO,GO,NT>(*B, exports, numPacketsPerLID(),
          exportLIDs(), constantNumPackets, distor, Tpetra::REPLACE);
      local_op_ok = 0;
    } catch (std::exception& e) {
      // This method should throw because it is not finished!
      local_op_ok = 1;
    }

    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsGraphAndCombine reported an error!" << endl;
      }
      gathervPrint(out, msg.str(), *comm);
      return; // no point in continuing
    }
  }

  // The test below uses the host Tpetra::CrsGraph interface to
  // compare graph values.  Thus, we need to do a fence before
  // comparing graph values, in order to ensure that changes made on
  // device are visible on host.
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;
  execution_space().fence ();

  int lclNumErrors = 0;

  out << "Comparing graphs after unpackCrsGraphAndCombine "
    "with CombineMode=REPLACE" << endl;
  {
    std::ostringstream errStrm;
    for (LO lclRow=0; lclRow<num_loc_rows; ++lclRow) {
      ArrayView<const LO> A_indices;
      A->getLocalRowView(lclRow, A_indices);

      ArrayView<const LO> B_indices;
      B->getLocalRowView(lclRow, B_indices);

      continue;
      /*
       * Test to be uncommented when unpackCrsGraphAndCombine is finished.
       *
      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        if (A_indices[i] != B_indices[i]) {
          errStrm << "ERROR: Proc " << world_rank << ", row " << lclRow
                  << ", A[" << i << "]=" << A_indices[i] << ", but "
                  <<   "B[" << i << "]=" << B_indices[i] << "!\n";
          ++curNumErrors;
        }
      }
      lclNumErrors += curNumErrors;
      */
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
    TEST_EQUALITY_CONST( gblNumErrors, 0 );
    if (gblNumErrors != 0) {
      if (world_rank == 0) {
        out << "unpackCrsGraphAndCombine comparison found " << gblNumErrors
            << " error" << (gblNumErrors != 1 ? "s" : "") << "!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    lclSuccess = success ? 1 : 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsGraphAndCombine comparison claims zero errors, "
          "but success is false on at least one process!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }
#endif // 0
}

// PackWithError sends intentionally bad inputs to pack/unpack to make sure
// that CrsGraph will detect the bad inputs and return the correct
// error diagnostics.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, PackWithError, LO, GO, NT)
{

  typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
  typedef typename crs_graph_type::packet_type packet_type;

  RCP<const Comm<int> > comm = getDefaultComm();
  const int world_rank = comm->getRank();

  out << "Creating graph" << endl;

  auto A = generate_test_graph<crs_graph_type>(comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Calling packCrsGraph" << endl;

  // Prepare arguments for pack.  This test is similar to the
  // PackThenUnpackAndCombine test,
  // but incorrect exportLIDs are sent in to induce a packing error.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for error string local to each process

  {
    LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
    Array<LO> exportLIDs(num_loc_rows);
    // exportLIDs[i] should equal i, but we set it to i+2
    for (LO i=0; i<num_loc_rows; i++) {
      exportLIDs[i] = i + 2;
    }

    Array<packet_type> exports;
    Array<size_t> numPacketsPerLID(num_loc_rows, 0);
    size_t constantNumPackets;
    {
      int local_op_ok;
      std::ostringstream msg;
      try {
        packCrsGraph<LO,GO,NT>(*A, exports, numPacketsPerLID(), exportLIDs(),
            constantNumPackets);
        local_op_ok = 1;
      } catch (std::exception& e) {
        local_op_ok = 0;
        msg << e.what();
      }
      if (local_op_ok == 1) {
        // Local pack should not be OK!  We requested bad local IDs be exported!
        errStrm << "Proc " << world_rank
                << ": packCrsGraph returned OK, but bad local IDs were requested!"
                << endl;
        lclSuccess = 0;
      }
    }
  }

  Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "packCrsGraph failed to notice bad export IDs on some process!" << endl;
    gathervPrint (out, errStrm.str (), *comm);
  }

  {
    // Let's try this again, but send in the wrong number of exportLIDs
    LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
    // Note the -1!
    out << "Allocating ids... ";
    Array<LO> exportLIDs(num_loc_rows-1);
    for (LO i=0; i < num_loc_rows-1; ++i) {
      exportLIDs[i] = i;
    }
    out << "done" << endl;

    Array<packet_type> exports;
    Array<size_t> numPacketsPerLID(num_loc_rows, 0);
    size_t constantNumPackets;
    out << "Calling packCrsGraph" << endl;
    {
      int local_op_ok;
      std::ostringstream msg;
      try {
        packCrsGraph<LO,GO,NT>(*A, exports, numPacketsPerLID(), exportLIDs(),
            constantNumPackets);
        local_op_ok = 1;
      } catch (std::exception& e) {
        local_op_ok = 0;
        msg << e.what();
      }
      if (local_op_ok == 1) {
        // Local pack should not be OK!  We requested too few local IDs be exported!
        errStrm << "Proc " << world_rank
                << ": packCrsGraph returned OK, but too few IDs "
                << "were requested to be exported!"
                << endl;
        lclSuccess = 0;
      }
    }
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }
}

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, PackThenUnpackAndCombine, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, PackWithError, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

} // namespace (anonymous)
