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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <map>
#include <random>
#include <vector>
#include <algorithm>
#include "Kokkos_DefaultNode.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import_Util2.hpp"

namespace {

template<class DeviceViewType>
struct create_views {
  DeviceViewType d;
  typename DeviceViewType::HostMirror h;
  create_views(const std::vector<typename DeviceViewType::value_type>& x)
  {
    using Kokkos::RangePolicy;
    using Kokkos::HostSpace;
    typedef RangePolicy<HostSpace::execution_space> range_policy;
    d = DeviceViewType("x_d", x.size());
    h = Kokkos::create_mirror_view(d);
    Kokkos::parallel_for(range_policy(0, x.size()),
        KOKKOS_LAMBDA (const int &i) {h(i) = x[i];});
    copy_to_device();
  }
  void copy_to_host() { Kokkos::deep_copy(h, d); }
  void copy_to_device() { Kokkos::deep_copy(d, h); }
};


template<class Ordinal, class Scalar>
void
generate_crs_entries(std::vector<size_t>& rowptr,
                     std::vector<size_t>& rowptr2,
                     std::vector<Ordinal>& colind,
                     std::vector<Ordinal>& colind2,
                     std::vector<Scalar>& vals,
                     std::vector<Scalar>& vals2,
                     int max_num_entries_per_row,
                     int num_cols)
{
  typedef typename std::vector<Ordinal>::size_type size_type;

  TEUCHOS_TEST_FOR_EXCEPTION(max_num_entries_per_row % 2 == 0,
    std::logic_error, "max_num_entries_per_row must be an odd integer!");

  TEUCHOS_TEST_FOR_EXCEPTION(num_cols % 2 == 0,
    std::logic_error, "num_cols must be an odd integer!");

  // Fill the CRS arrays, use random values
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<Scalar> dist(1.0, 2.0);
  int row = 0;
  while (true) {
    int m = (max_num_entries_per_row - 1) / 2;
    int start = std::max(0, row - m);
    int num_cols_this_row = std::min(row + 1 + m,
        std::min(max_num_entries_per_row, m+num_cols-row));
    int end = start + num_cols_this_row;

    rowptr.push_back(static_cast<size_t>(colind.size()));
    rowptr2.push_back(static_cast<size_t>(colind2.size()));

    for (int col=start; col<end; col++) {
      colind.push_back(col);
      colind2.push_back(col);

      Scalar rands = dist(gen);
      vals.push_back(rands);
      vals2.push_back(rands);

      // Create duplicate for merge test
      colind2.push_back(col);
      vals2.push_back(rands);
    }

    if (row > 1 && num_cols_this_row == 1 + m) break;
    row++;

  }
  size_type num_entries = colind.size();
  rowptr.push_back(static_cast<size_t>(num_entries));

  size_type num_entries2 = colind2.size();
  rowptr2.push_back(static_cast<size_t>(num_entries2));

  return;
}


template<class Ordinal, class Scalar>
void
shuffle_crs_entries(std::vector<Ordinal>& colind_rand,
                    std::vector<Scalar>& vals_rand,
                    const std::vector<size_t>& rowptr,
                    const std::vector<Ordinal>& colind,
                    const std::vector<Scalar>& vals)
{

  typedef typename std::vector<Ordinal> colind_type;

  // Randomly shuffle values and column indices
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<Scalar> dist(1.0, 2.0);

  for (size_t i=0; i < static_cast<size_t>(rowptr.size()-1); ++i) {
    size_t num_cols_this_row = rowptr[i+1] - rowptr[i];
    colind_type ix(num_cols_this_row);
    std::iota(ix.begin(), ix.end(), rowptr[i]);
    std::shuffle(ix.begin(), ix.end(), gen);
    for (size_t j=rowptr[i], k=0; j<rowptr[i+1]; ++j, ++k) {
      vals_rand[j] = vals[ix[k]];
      colind_rand[j] = colind[ix[k]];
    }
  }

  return;
}


//
// UNIT TESTS
//


// Unit Test the functionality in Tpetra_Import_Util
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Import_Util, SortCrsEntries, Scalar, LO, GO)
{

  typedef typename std::vector<size_t>  rowptr_type;
  typedef typename std::vector<GO>  colind_type;
  typedef typename std::vector<Scalar> vals_type;
  typedef typename colind_type::size_type size_type;

  int max_num_entries_per_row = 7;  // should be odd
  int num_cols = 15; // should be odd

  rowptr_type  rowptr, rowptr2;
  colind_type  colind, colind2;
  vals_type vals, vals2;

  generate_crs_entries(rowptr, rowptr2, colind, colind2, vals, vals2,
                       max_num_entries_per_row, num_cols);

  {

    //
    // Sort the CRS entries by column index
    //

    size_type num_entries = colind.size();

    // Randomly shuffle values and column indices
    vals_type  vals_rand(num_entries);
    colind_type  colind_rand(num_entries);
    shuffle_crs_entries(colind_rand, vals_rand, rowptr, colind, vals);

    // Make copy for sortCrsEntries w/o values
    colind_type colind_rand_copy(colind_rand.begin(), colind_rand.end());

    //
    // Sort the GIDs and associated values
    //

    // Create array views of data
    auto rowptr_av = Teuchos::ArrayView<size_t>(rowptr);
    auto colind_rand_av = Teuchos::ArrayView<GO>(colind_rand);
    auto vals_rand_av = Teuchos::ArrayView<Scalar>(vals_rand);

    // Do the sort
    Tpetra::Import_Util::sortCrsEntries<Scalar, GO>(rowptr_av, colind_rand_av, vals_rand_av);

    // At this point, colind_rand and vals_rand should be sorted and equal to
    // their original arrays, respectively.
    TEST_COMPARE_ARRAYS(colind, colind_rand);
    TEST_COMPARE_FLOATING_ARRAYS(vals, vals_rand, 1.e-12);

    //
    // Sort the GIDs w/o values
    //

    // Create array views of data
    auto colind_rand_copy_av = Teuchos::ArrayView<GO>(colind_rand_copy);

    // Do the sort
    Tpetra::Import_Util::sortCrsEntries<GO>(rowptr_av, colind_rand_copy_av);

    // At this point, colind_rand_copy should be sorted and equal to their
    // original array
    TEST_COMPARE_ARRAYS(colind, colind_rand_copy);

  }

  {
    //
    // Sort and merge the CRS entries by column index
    //

    size_type num_entries = colind2.size();

    // Randomly shuffle values and column indices for merge test
    vals_type  vals_rand(num_entries);
    colind_type  colind_rand(num_entries);
    shuffle_crs_entries(colind_rand, vals_rand, rowptr2, colind2, vals2);

    // Make copies for sortAndMergeCrsEntries w/o values
    rowptr_type rowptr2_copy(rowptr2.begin(), rowptr2.end());
    colind_type colind_rand_copy(colind_rand.begin(), colind_rand.end());

    //
    // Sort the GIDs and associated values
    //

    // Create array views of data
    auto rowptr_av = Teuchos::ArrayView<size_t>(rowptr2);
    auto colind_rand_av = Teuchos::ArrayView<GO>(colind_rand);
    auto vals_rand_av = Teuchos::ArrayView<Scalar>(vals_rand);

    // Do the sort
    Tpetra::Import_Util::sortAndMergeCrsEntries<Scalar, GO>(rowptr_av, colind_rand_av, vals_rand_av);

    // At this point, the row pointers and column indices should be the same as
    // the version above due to the merge/shrink
    size_type new_num_entries = rowptr2[rowptr2.size()-1];
    TEST_EQUALITY(colind.size(), new_num_entries);

    TEST_COMPARE_ARRAYS(rowptr, rowptr2);

    colind_rand.resize(new_num_entries);
    TEST_COMPARE_ARRAYS(colind, colind_rand);

    // At this point, the values should be twice the vals above due to the
    // merge/shrink
    vals2.resize(new_num_entries);
    vals_rand.resize(new_num_entries);
    for (size_type i=0; i<new_num_entries; i++) vals2[i] = 2.*vals[i];
    TEST_COMPARE_FLOATING_ARRAYS(vals2, vals_rand, 1.e-12);

    //
    // Sort the GIDs w/o values
    //

    // Create array views of data
    auto rowptr_copy_av = Teuchos::ArrayView<size_t>(rowptr2_copy);
    auto colind_rand_copy_av = Teuchos::ArrayView<GO>(colind_rand_copy);

    // Do the sort
    Tpetra::Import_Util::sortAndMergeCrsEntries<GO>(rowptr_copy_av, colind_rand_copy_av);

    // At this point, the row pointers and column indices should be the same as
    // the version above due to the merge/shrink
    new_num_entries = rowptr2_copy[rowptr2_copy.size()-1];
    TEST_EQUALITY(colind.size(), new_num_entries);

    TEST_COMPARE_ARRAYS(rowptr, rowptr2_copy);

    colind_rand_copy.resize(new_num_entries);
    TEST_COMPARE_ARRAYS(colind, colind_rand_copy);

  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Import_Util, SortCrsEntriesKokkos, Scalar, LO, GO, NT)
{

  typedef typename Tpetra::CrsMatrix<Scalar,LO,GO,NT>::local_matrix_type local_matrix_type;
  typedef typename local_matrix_type::StaticCrsGraphType graph_type;
  typedef typename graph_type::row_map_type::non_const_type rowptr_view_type;
  typedef typename graph_type::entries_type::non_const_type colind_view_type;
  typedef typename local_matrix_type::values_type::non_const_type vals_view_type;

  typedef typename std::vector<typename rowptr_view_type::value_type> rowptr_type;
  typedef typename std::vector<typename colind_view_type::value_type> colind_type;
  typedef typename std::vector<typename vals_view_type::value_type> vals_type;
  typedef typename colind_type::size_type size_type;

  // Map is not actually used, but is needed to instantiate Kokkos
  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const Tpetra::global_size_t INVALID =
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  auto dummy_map = Tpetra::Map<LO,GO,NT>(INVALID, 1, 0, comm);

  int max_num_entries_per_row = 7;  // should be odd
  int num_cols = 15; // should be odd

  rowptr_type  rowptr, rowptr2;
  colind_type  colind, colind2;
  vals_type vals, vals2;

  generate_crs_entries(rowptr, rowptr2, colind, colind2, vals, vals2,
                       max_num_entries_per_row, num_cols);

  {
    //
    // Sort the CRS entries by column index
    //

    size_type num_entries = colind.size();

    // Randomly shuffle values and column indices
    vals_type  vals_rand(num_entries);
    colind_type  colind_rand(num_entries);
    shuffle_crs_entries(colind_rand, vals_rand, rowptr, colind, vals);

    // Create mirror views of the CRS entries
    auto rowptr_views = create_views<rowptr_view_type>(rowptr);
    auto colind_rand_views = create_views<colind_view_type>(colind_rand);
    auto colind_rand_copy_views = create_views<colind_view_type>(colind_rand);
    auto vals_rand_views = create_views<vals_view_type>(vals_rand);

    //
    // Sort the GIDs and associated values
    //
    Tpetra::Import_Util::sortCrsEntries<rowptr_view_type, colind_view_type,
      vals_view_type>(rowptr_views.d, colind_rand_views.d, vals_rand_views.d);

    // Copy back to host
    colind_rand_views.copy_to_host();
    vals_rand_views.copy_to_host();

    // At this point, colind_rand and vals_rand should be sorted and equal to
    // their original arrays, respectively.
    TEST_COMPARE_ARRAYS(colind, colind_rand_views.h);
    TEST_COMPARE_FLOATING_ARRAYS(vals, vals_rand_views.h, 1.e-12);

    //
    // Sort the GIDs w/o values
    //

    // Do the sort
    Tpetra::Import_Util::sortCrsEntries<rowptr_view_type, colind_view_type>(
        rowptr_views.d, colind_rand_copy_views.d);

    // Copy back to host
    colind_rand_copy_views.copy_to_host();

    // At this point, colind_rand_copy should be sorted and equal to their
    // original array
    TEST_COMPARE_ARRAYS(colind, colind_rand_copy_views.h);

  }

}


  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Import_Util, SortCrsEntries, SC, LO, GO )

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Import_Util, SortCrsEntriesKokkos, SC, LO, GO, NT )

  // Note: This test fails.  Should fix later.
  //      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ReverseImportExport, doImport, ORDINAL, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Test CrsMatrix for all Scalar, LO, GO template parameter
  // combinations, and the default Node type.
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO )

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )

}
