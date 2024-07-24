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
#include "Tpetra_Details_crsUtils.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <random>
#include <set>

namespace { // anonymous


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


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, PackThenUnpackAndCombine, LO, GO, NT)
{

  // This is a relatively simple test.  We wish to create a tridiagonal graph
  // A by the following process:
  //
  // 1. Create diagonal graph A
  // 2. Create graph B containing only off-diagonal values
  // 3. Pack B and unpack entries in to A
  //
  // The end result is the tridiagonal graph we originally wanted.
  //
  // This test exercises CrsGraph's ability to unpack in to a graph with static
  // profile.
  using Teuchos::tuple;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using packet_type = typename graph_type::packet_type;
//  using import_type = Tpetra::Import<LO,GO,NT>;

  auto comm = getDefaultComm();
  const auto my_rank = comm->getRank();

  // Create diagonal graph.
  auto num_loc_rows = static_cast<LO>(4);
  const auto num_gbl_rows = Tpetra::global_size_t(num_loc_rows*comm->getSize());
  auto map1 = rcp(new map_type(num_gbl_rows, 0, comm));
  auto A = rcp(new graph_type(map1, 1));
  for (LO loc_row=0; loc_row<num_loc_rows; loc_row++) {
    const auto gbl_row = map1->getGlobalElement(loc_row);
    A->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row));
  }

  // Off diagonal graph with half-bandwidth=1 and no diagonal entries
  out << "Building second graph" << endl;
  auto B = rcp(new graph_type(map1, 2)); 
  for (LO loc_row=0; loc_row<num_loc_rows; loc_row++) {
    const auto gbl_row = map1->getGlobalElement(loc_row);
    // B[0,0:1] = [-, 1]
    if (gbl_row == 0)
      B->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row+1));
    // B[N-1,N-2:N-1] = [1, -]
    else if (static_cast<Tpetra::global_size_t>(gbl_row) == num_gbl_rows-1)
      B->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row-1));
    // B[I,I-1:I+1] = [1, -, 1]
    else
      B->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row-1, gbl_row+1));
  }
  B->fillComplete();

  auto loc_success = 1; // to be revised below
  auto gbl_success = 0; // output argument

  out << "Preparing arguments for packCrsGraph" << endl;

  auto export_lids = Array<LO>(num_loc_rows); // input argument
  for (auto i=0; i<num_loc_rows; ++i) {
    export_lids[i] = static_cast<LO>(i); // pack all the rows
  }
  auto exports = Array<packet_type>(); // output argument; to be realloc'd
  auto num_packets_per_lid = Array<size_t>(num_loc_rows, 0); // output argument
  size_t const_num_packets; // output argument

  out << "Calling packCrsGraph" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsGraph<LO,GO,NT>(*B, exports, num_packets_per_lid(), export_lids(),
          const_num_packets);
      local_op_ok = 1;
    } catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }
    TEST_ASSERT(local_op_ok == 1);
    loc_success = success ? 1 : 0;
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MIN, loc_success, outArg(gbl_success));
    TEST_EQUALITY_CONST(gbl_success, 1);
    if (gbl_success != 1) {
      if (my_rank == 0) {
        out << "packCrsGraph reported an error!" << endl;
      }
      gathervPrint(out, msg.str(), *comm);
      out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  return;
#if 0

  // Now unpack in to the static graph
  out << "Calling unpackCrsGraphAndCombine" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    unpackCrsGraphAndCombine<LO,GO,NT>(*A, exports, num_packets_per_lid(),
        export_lids(), const_num_packets, Tpetra::REPLACE);
    local_op_ok = 1;

    TEST_ASSERT(local_op_ok == 1);
    loc_success = success ? 1 : 0;
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MIN, loc_success, outArg(gbl_success));
    TEST_EQUALITY_CONST(gbl_success, 1);
    if (gbl_success != 1) {
      if (my_rank == 0) {
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
  A->fillComplete();
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;
  using gids_type = typename graph_type::nonconst_global_inds_host_view_type;
  execution_space().fence ();

  auto loc_num_errs = 0;

  out << "Comparing graphs after unpackCrsGraphAndCombine" << endl;
  {
    std::ostringstream errStrm;
    for (LO loc_row=0; loc_row<num_loc_rows; ++loc_row) {
      const auto gbl_row = map1->getGlobalElement(loc_row);
      size_t num_entries = 3;
      gids_type A_indices(num_entries);
      A->getGlobalRowCopy(gbl_row, A_indices, num_entries);
      Tpetra::sort(A_indices, num_entries);

      auto errors = 0; // Herb Sutter loves you :)
      if (gbl_row == 0) {
        // A[0,0:1] = [1, 1]
        if (num_entries != 2) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", expected row to have 2 indices not " << num_entries << "!\n";
          ++errors;
        } else if (!(A_indices[0]==gbl_row && A_indices[1]==gbl_row+1)) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", incorrect indices: "
                  << A_indices(0,num_entries)
                  << " != [" << gbl_row << ", " << gbl_row+1 << "]\n";
          ++errors;
        }
      } else if (static_cast<Tpetra::global_size_t>(gbl_row) == num_gbl_rows-1) {
        // B[N-1,N-2:N-1] = [1, -]
        if (num_entries != 2) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", expected row to have 3 indices not " << num_entries << "!\n";
          ++errors;
        } else if (!(A_indices[0]==gbl_row-1 && A_indices[1]==gbl_row)) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", incorrect indices: "
                  << A_indices(0,num_entries)
                  << " != [" << gbl_row-1 << ", " << gbl_row << "]\n";
          ++errors;
        }
      }
      else {
        // B[I,I-1:I+1] = [1, -, 1]
        if (num_entries != 3) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", expected row to have 3 indices not " << num_entries << "!\n";
          ++errors;
        } else if (!(A_indices[0]==gbl_row-1 && A_indices[1]==gbl_row && A_indices[2]==gbl_row+1)) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << gbl_row
                  << ", incorrect indices: "
                  << A_indices(0,num_entries)
                  << " != [" << gbl_row-1 << ", " << gbl_row << ", " << gbl_row+1 << "]\n";
          ++errors;
        }
      }
      loc_num_errs += errors;
    }
    TEST_ASSERT(loc_num_errs == 0);

    auto gbl_num_errs = 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_SUM, loc_num_errs, outArg(gbl_num_errs));
    TEST_EQUALITY_CONST(gbl_num_errs, 0);
    if (gbl_num_errs != 0) {
      if (my_rank == 0) {
        out << "unpackCrsGraphAndCombine comparison found " << gbl_num_errs
            << " error" << (gbl_num_errs != 1 ? "s" : "") << "!" << endl;
      }
      gathervPrint(out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    loc_success = success ? 1 : 0;
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MIN, loc_success, outArg(gbl_success));
    TEST_EQUALITY_CONST(gbl_success, 1 );
    if (gbl_success != 1) {
      if (my_rank == 0) {
        out << "unpackCrsGraphAndCombine comparison claims zero errors, "
          "but success is false on at least one process!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }
#endif // 0
}


#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, PackThenUnpackAndCombine, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
