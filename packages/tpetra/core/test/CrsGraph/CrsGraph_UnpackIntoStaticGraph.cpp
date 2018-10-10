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
#include "Tpetra_Details_resizeRowPtrs.hpp"
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
using Tpetra::Details::resizeRowPtrsAndIndices;
using std::endl;

template<class T>
bool
essentially_equal(T a, T b) {
  typedef Kokkos::ArithTraits<T> KAT;
  const auto eps = KAT::eps();
  return KAT::abs(a - b) <= ( (KAT::abs(a) > KAT::abs(b) ? KAT::abs(b) : KAT::abs(a)) * eps);
}


TEUCHOS_UNIT_TEST(CrsGraph, ResizeRowPointersAndIndices)
{

  int argc = 0;
  char* argv[] = {NULL};
  char** argv2 = &argv[0];
  Tpetra::ScopeGuard tpetra_scope(&argc, &argv2);

  using ordinal_type = size_t;
  using view_type = Kokkos::View<ordinal_type*>;

  ordinal_type num_row = 4;
  ordinal_type num_indices_per_row = 5;
  ordinal_type num_indices = num_indices_per_row * num_row;
  auto row_ptrs_beg = view_type("beg", num_row+1);
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

  auto num_packets_per_lid = view_type("num packets", num_row);
  for (ordinal_type i=0; i<num_row; i++) num_packets_per_lid(i) = i;
  ordinal_type num_extra = num_row*(num_packets_per_lid(0) + num_packets_per_lid(num_row-1))/2;

  resizeRowPtrsAndIndices<view_type, view_type, view_type>(row_ptrs_beg, row_ptrs_end, indices,
                                                           num_packets_per_lid, false);

  TEST_ASSERT(indices.size() == num_indices + num_extra);

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
  auto D = rcp(new graph_type(map1, 1, Tpetra::StaticProfile));
  for (LO loc_row=0; loc_row<num_loc_rows; loc_row++) {
    const auto gbl_row = map1->getGlobalElement(loc_row);
    D->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row));
  }
  D->fillComplete();

  // Off diagonal graph with bandwidth=1
  auto X = rcp(new graph_type(map1, 2));
  for (LO loc_row=0; loc_row<num_loc_rows; loc_row++) {
    const auto gbl_row = map1->getGlobalElement(loc_row);
    // X[0,0:1] = [-, 1]
    if (gbl_row == 0)
      X->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row));
    // X[N-1,N-2:N-1] = [1, -]
    else if (static_cast<Tpetra::global_size_t>(gbl_row) == num_gbl_rows-1)
      X->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row-1));
    // X[I,I-1:I+1] = [1, -, 1]
    else
      X->insertGlobalIndices(gbl_row, tuple<GO>(gbl_row-1, gbl_row+1));
  }
  X->fillComplete();

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
  auto distor = Tpetra::Distributor(comm); // argument required, but not used

  out << "Calling packCrsGraph" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsGraph<LO,GO,NT>(*X, exports, num_packets_per_lid(), export_lids(),
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

  // Now make sure that the pack is correct by creating an empty graph and
  // unpacking in to it.  The graph should end up being the same as the above graph.
  out << "Building second graph" << endl;

#ifdef KOKKOS_ENABLE_SERIAL
  const bool atomic_updates = ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
  const bool atomic_updates = true;
#endif // KOKKOS_ENABLE_SERIAL

  out << "Calling unpackCrsGraphAndCombine with "
      << "CombineMode=Tpetra::REPLACE" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    D->resumeFill();
    unpackCrsGraphAndCombine<LO,GO,NT>(*D, exports, num_packets_per_lid(),
        export_lids(), const_num_packets, distor, Tpetra::REPLACE, atomic_updates);
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
  execution_space::fence ();

  auto loc_num_errs = 0;

  out << "Comparing graphs after unpackCrsGraphAndCombine "
    "with CombineMode=REPLACE" << endl;
  {
    std::ostringstream errStrm;
    for (LO loc_row=0; loc_row<num_loc_rows; ++loc_row) {
      ArrayView<const LO> A_indices;
      D->getLocalRowView(loc_row, A_indices);

      //ArrayView<const LO> B_indices;
      //B->getLocalRowView(loc_row, B_indices);

      continue;
      /*
       * Test to be uncommented when unpackCrsGraphAndCombine is finished.
       *
      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        if (A_indices[i] != B_indices[i]) {
          errStrm << "ERROR: Proc " << my_rank << ", row " << loc_row
                  << ", A[" << i << "]=" << A_indices[i] << ", but "
                  <<   "B[" << i << "]=" << B_indices[i] << "!\n";
          ++curNumErrors;
        }
      }
      loc_num_errs += curNumErrors;
      */
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
