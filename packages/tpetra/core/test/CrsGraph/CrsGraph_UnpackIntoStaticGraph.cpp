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

#include "Tpetra_TestingUtilities.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_packCrsGraph.hpp"
#include "Tpetra_Details_unpackCrsGraphAndCombine.hpp"
#include "Tpetra_Details_padCrsArrays.hpp"
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
using Tpetra::Details::unpackCrsGraphAndCombine;
using Tpetra::Details::padCrsArrays;
using std::endl;

TEUCHOS_UNIT_TEST(CrsGraph, ResizeRowPointersAndIndices_1)
{
  using device_type = typename Tpetra::Map<>::device_type;
  using execution_space = typename device_type::execution_space;
  using ordinal_type = size_t;
  using view_type = Kokkos::View<ordinal_type*, device_type>;
  using size_type = typename view_type::size_type;

  ordinal_type num_row = 4;
  ordinal_type num_indices_per_row = 5;
  ordinal_type num_indices = num_indices_per_row * num_row;
  auto row_ptrs_beg = view_type("beg", num_row+1);
  // this assumes UVM
  for (ordinal_type i=0; i<num_row+1; i++) row_ptrs_beg(i) = num_indices_per_row*i;

  auto row_ptrs_end = view_type("end", num_row);
  for (ordinal_type i=0; i<num_row; i++) row_ptrs_end(i) = row_ptrs_beg(i+1);

  auto indices = view_type("indices", num_indices);
  for (ordinal_type i=0; i<num_row; i++) {
    auto start = row_ptrs_beg(i);
    auto end = row_ptrs_beg(i+1);
    for (ordinal_type j=start, k=0; j<end; j++, k++) {
      indices(j) = (i + 1) * num_indices_per_row + k;
    }
  }

  auto import_lids = view_type("import lids", num_row);
  auto num_packets_per_lid = view_type("num packets", num_row);
  for (ordinal_type i=0; i<num_row; i++) {
   import_lids(i) = i;
   num_packets_per_lid(i) = i;
  }
  ordinal_type num_extra = num_row*(num_packets_per_lid(0) + num_packets_per_lid(num_row-1))/2;

  Kokkos::UnorderedMap<ordinal_type,ordinal_type,device_type> padding(import_lids.size());
  execution_space::fence();
  for (size_type i=0; i<import_lids.size(); i++){
    padding.insert(import_lids(i), num_packets_per_lid(i));
  }
  execution_space::fence();
  TEST_ASSERT(!padding.failed_insert());

  padCrsArrays(row_ptrs_beg, row_ptrs_end, indices, padding);
  TEST_ASSERT(indices.size() == static_cast<size_type>(num_indices + num_extra));

  {
    // make sure row_ptrs_beg is where it should be
    bool offsets_ok = row_ptrs_beg(0) == 0;
    for (ordinal_type i=1; i<num_row; i++) {
      auto expected = row_ptrs_beg(i-1) + num_indices_per_row + num_packets_per_lid(i-1);
      if (row_ptrs_beg(i) != expected) offsets_ok = false;
    }
    if (offsets_ok) offsets_ok = row_ptrs_beg(num_row) == indices.size();
    TEST_ASSERT(offsets_ok);
  }

  {
    // make sure indices were shifted correctly
    bool indices_ok = true;
    for (ordinal_type i=0; i<num_row; i++) {
      auto start = row_ptrs_beg(i);
      auto end = row_ptrs_end(i);
      for (ordinal_type j=start, k=0; j<end; j++, k++) {
        auto expected = (i + 1) * num_indices_per_row + k;
        if (expected != indices(j)) indices_ok = false;
      }
    }
    TEST_ASSERT(indices_ok);
  }
}

TEUCHOS_UNIT_TEST(CrsGraph, ResizeRowPointersAndIndices_2)
{
  typedef typename Tpetra::Map<>::device_type device_type;
  using execution_space = typename device_type::execution_space;
  using ordinal_type = size_t;
  using view_type = Kokkos::View<ordinal_type*, device_type>;
  using size_type = typename view_type::size_type;

  auto row_ptrs_beg = view_type("beg", 4);
  auto row_ptrs_end = view_type("end", 3);
  auto indices = view_type("indices", 9);

  // Row 1, 3 allocations, 3 used
  row_ptrs_beg(0) = 0; row_ptrs_end(0) = 3;
  indices(0) = 1; indices(1) = 2; indices(2) = 3;  // Row 1

  // Row 2, 3 allocations only 1 used
  row_ptrs_beg(1) = 3; row_ptrs_end(1) = 4;
  indices(3) = 4;

  // Row 3, 3 allocations, only 2 used
  row_ptrs_beg(2) = 6; row_ptrs_end(2) = 8;
  indices(6) = 7; indices(7) = 8;

  row_ptrs_beg(3) = 9;

  // Import 5 extra values in to Row 1 and 3 extra in to Row 3
  auto import_lids = view_type("import lids", 2);
  auto num_packets_per_lid = view_type("num packets", 2);

  // Import LIDs not ordered
  import_lids(0) = 2;
  num_packets_per_lid(0) = 3;

  import_lids(1) = 0;
  num_packets_per_lid(1) = 5;

  Kokkos::UnorderedMap<ordinal_type,ordinal_type,device_type> padding(import_lids.size());
  execution_space::fence();
  for (size_type i=0; i<import_lids.size(); i++){
    padding.insert(import_lids(i), num_packets_per_lid(i));
  }
  execution_space::fence();
  TEST_ASSERT(!padding.failed_insert());
  padCrsArrays(row_ptrs_beg, row_ptrs_end, indices, padding);

  // Check row offsets
  TEST_ASSERT(row_ptrs_beg(0) == 0);
  TEST_ASSERT(row_ptrs_beg(1) == 8);
  TEST_ASSERT(row_ptrs_beg(2) == 11);
  TEST_ASSERT(row_ptrs_beg(3) == 16);
  TEST_ASSERT(indices.size() == 16);

  TEST_ASSERT(row_ptrs_end(0) == 3);
  TEST_ASSERT(row_ptrs_end(1) == 9);
  TEST_ASSERT(row_ptrs_end(2) == 13);

  // Row 1
  TEST_ASSERT(indices(0) == 1);
  TEST_ASSERT(indices(1) == 2);
  TEST_ASSERT(indices(2) == 3);
  // 5 extra entries in Row 1

  // Row 2
  TEST_ASSERT(indices(8) == 4);
  // 2 unused indices in Row 2

  // Row 2
  TEST_ASSERT(indices(11) == 7);
  TEST_ASSERT(indices(12) == 8);

}

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
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;
//  using import_type = Tpetra::Import<LO,GO,NT>;

  auto comm = getDefaultComm();
  const auto my_rank = comm->getRank();

  // Create diagonal graph.
  auto num_loc_rows = static_cast<LO>(4);
  const auto num_gbl_rows = Tpetra::global_size_t(num_loc_rows*comm->getSize());
  auto map1 = rcp(new map_type(num_gbl_rows, 0, comm));
  auto A = rcp(new graph_type(map1, 1, Tpetra::StaticProfile));
  for (LO loc_row=0; loc_row<num_loc_rows; loc_row++) {
    const auto gbl_row = map1->getGlobalElement(loc_row);
    A->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row));
  }

  // Off diagonal graph with half-bandwidth=1 and no diagonal entries
  out << "Building second graph" << endl;
  auto B = rcp(new graph_type(map1, 2)); // could use StaticProfile
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

  // We're not actually communicating in this test; we just need the Distributor
  // for the interface of packCrsGraph (which doesn't use it).  Consider changing
  // packCrsGraph's interface so it doesn't take a Distributor?  No, because
  // Distributor has index permutation information that we could use to pack in
  // a particular order and thus avoid the slow path in Distributor::doPosts.
  auto distor = Tpetra::Distributor(comm);

  out << "Calling packCrsGraph" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsGraph<LO,GO,NT>(*B, exports, num_packets_per_lid(), export_lids(),
          const_num_packets, distor);
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

  // Now unpack in to the static graph
  out << "Calling unpackCrsGraphAndCombine" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    unpackCrsGraphAndCombine<LO,GO,NT>(*A, exports, num_packets_per_lid(),
        export_lids(), const_num_packets, distor, Tpetra::REPLACE);
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
  execution_space::fence ();

  auto loc_num_errs = 0;

  out << "Comparing graphs after unpackCrsGraphAndCombine" << endl;
  {
    std::ostringstream errStrm;
    for (LO loc_row=0; loc_row<num_loc_rows; ++loc_row) {
      const auto gbl_row = map1->getGlobalElement(loc_row);
      size_t num_entries = 3;
      Array<GO> A_indices(num_entries);
      A->getGlobalRowCopy(gbl_row, A_indices(), num_entries);
      std::sort(A_indices.begin(), A_indices.begin()+num_entries);

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
