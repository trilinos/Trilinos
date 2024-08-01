// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <map>
#include <random>
#include <vector>
#include <algorithm>
#include "Tpetra_KokkosCompat_DefaultNode.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace {

template<class DeviceViewType>
struct create_device_views {
  DeviceViewType d;
  typename DeviceViewType::HostMirror h;
  create_device_views(
    const std::vector<typename DeviceViewType::value_type>& x
  )
  {
    d = DeviceViewType("x_d", x.size());
    h = Kokkos::create_mirror_view(d);
    for (size_t i=0; i<static_cast<size_t>(x.size()); i++) h(i) = x[i];
    copy_to_device();
  }
  void copy_to_host() { Kokkos::deep_copy(h, d); }
  void copy_to_device() { Kokkos::deep_copy(d, h); }
};


template<class scalar_type, class ordinal_type, class index_type>
void
generate_crs_entries(std::vector<index_type>& rowptr,
                     std::vector<index_type>& rowptr2,
                     std::vector<ordinal_type>& colind,
                     std::vector<ordinal_type>& colind2,
                     std::vector<scalar_type>& vals,
                     std::vector<scalar_type>& vals2,
                     int max_num_entries_per_row,
                     int num_cols)
{
  typedef typename std::vector<ordinal_type>::size_type size_type;

  TEUCHOS_TEST_FOR_EXCEPTION(max_num_entries_per_row % 2 == 0,
    std::logic_error, "max_num_entries_per_row must be an odd integer!");

  TEUCHOS_TEST_FOR_EXCEPTION(num_cols % 2 == 0,
    std::logic_error, "num_cols must be an odd integer!");

  // Fill the CRS arrays, use random values
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<typename Kokkos::ArithTraits<scalar_type>::mag_type> dist(1.0, 2.0);
  int row = 0;
  while (true) {
    int m = (max_num_entries_per_row - 1) / 2;
    int start = std::max(0, row - m);
    int num_cols_this_row = std::min(row + 1 + m,
        std::min(max_num_entries_per_row, m+num_cols-row));
    int end = start + num_cols_this_row;

    rowptr.push_back(static_cast<index_type>(colind.size()));
    rowptr2.push_back(static_cast<index_type>(colind2.size()));

    for (int col=start; col<end; col++) {
      colind.push_back(col);
      colind2.push_back(col);

      scalar_type rands = static_cast<scalar_type> (dist(gen));
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
  rowptr.push_back(static_cast<index_type>(num_entries));

  size_type num_entries2 = colind2.size();
  rowptr2.push_back(static_cast<index_type>(num_entries2));

  return;
}


template<class scalar_type, class ordinal_type, class index_type>
void
shuffle_crs_entries(std::vector<ordinal_type>& colind_rand,
                    std::vector<scalar_type>& vals_rand,
                    const std::vector<index_type>& rowptr,
                    const std::vector<ordinal_type>& colind,
                    const std::vector<scalar_type>& vals)
{

  typedef typename std::vector<ordinal_type> colind_type;
  typedef typename std::vector<index_type>::size_type size_type;

  // Randomly shuffle values and column indices
  std::random_device rd;
  std::mt19937 gen(rd());

  for (size_type i=0; i < static_cast<size_type>(rowptr.size()-1); ++i) {
    size_type num_cols_this_row = rowptr[i+1] - rowptr[i];
    colind_type ix(num_cols_this_row);
    std::iota(ix.begin(), ix.end(), rowptr[i]);
    std::shuffle(ix.begin(), ix.end(), gen);
    for (size_type j=rowptr[i], k=0; j<rowptr[i+1]; ++j, ++k) {
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

  using Tpetra::Import_Util::sortCrsEntries;
  using Tpetra::Import_Util::sortAndMergeCrsEntries;

  typedef size_t index_type;
  typedef Scalar scalar_type;
  typedef GO ordinal_type;

  typedef typename std::vector<index_type>  rowptr_type;
  typedef typename std::vector<ordinal_type>  colind_type;
  typedef typename std::vector<scalar_type> vals_type;
  typedef typename colind_type::size_type size_type;

  int max_num_entries_per_row = 7;  // should be odd
  int num_cols = 15; // should be odd

  rowptr_type  rowptr, rowptr2;
  colind_type  colind, colind2;
  vals_type vals, vals2;

  generate_crs_entries<scalar_type,ordinal_type,index_type>(
      rowptr, rowptr2, colind, colind2, vals, vals2, max_num_entries_per_row, num_cols);

  {

    //
    // Sort the CRS entries by column index
    //

    size_type num_entries = colind.size();

    // Randomly shuffle values and column indices
    vals_type  vals_rand(num_entries);
    colind_type  colind_rand(num_entries);
    shuffle_crs_entries<scalar_type, ordinal_type, index_type>(
        colind_rand, vals_rand, rowptr, colind, vals);

    // Make copy for sortCrsEntries w/o values
    colind_type colind_rand_copy(colind_rand.begin(), colind_rand.end());

    //
    // Sort the GIDs and associated values
    //

    // Create array views of data
    auto rowptr_av = Teuchos::ArrayView<index_type>(rowptr);
    auto colind_rand_av = Teuchos::ArrayView<ordinal_type>(colind_rand);
    auto vals_rand_av = Teuchos::ArrayView<scalar_type>(vals_rand);

    // Do the sort
    sortCrsEntries<scalar_type, ordinal_type>(rowptr_av, colind_rand_av, vals_rand_av);

    // At this point, colind_rand and vals_rand should be sorted and equal to
    // their original arrays, respectively.
    TEST_COMPARE_ARRAYS(colind, colind_rand);
    TEST_COMPARE_FLOATING_ARRAYS(vals, vals_rand, 1.e-12);

    //
    // Sort the GIDs w/o values
    //

    // Create array views of data
    auto colind_rand_copy_av = Teuchos::ArrayView<ordinal_type>(colind_rand_copy);

    // Do the sort
    sortCrsEntries<ordinal_type>(rowptr_av, colind_rand_copy_av);

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
    shuffle_crs_entries<scalar_type,ordinal_type,index_type>(
        colind_rand, vals_rand, rowptr2, colind2, vals2);

    // Make copies for sortAndMergeCrsEntries w/o values
    rowptr_type rowptr2_copy(rowptr2.begin(), rowptr2.end());
    colind_type colind_rand_copy(colind_rand.begin(), colind_rand.end());

    //
    // Sort the GIDs and associated values
    //

    // Create array views of data
    auto rowptr_av = Teuchos::ArrayView<index_type>(rowptr2);
    auto colind_rand_av = Teuchos::ArrayView<ordinal_type>(colind_rand);
    auto vals_rand_av = Teuchos::ArrayView<scalar_type>(vals_rand);

    // Do the sort
    sortAndMergeCrsEntries<scalar_type,ordinal_type>(rowptr_av, colind_rand_av, vals_rand_av);

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
    for (size_type i=0; i<new_num_entries; i++) 
      vals2[i] = scalar_type(2.)*vals[i];
    TEST_COMPARE_FLOATING_ARRAYS(vals2, vals_rand, 1.e-12);

    //
    // Sort the GIDs w/o values
    //

    // Create array views of data
    auto rowptr_copy_av = Teuchos::ArrayView<index_type>(rowptr2_copy);
    auto colind_rand_copy_av = Teuchos::ArrayView<ordinal_type>(colind_rand_copy);

    // Do the sort
    sortAndMergeCrsEntries<ordinal_type>(rowptr_copy_av, colind_rand_copy_av);

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
  using Tpetra::Import_Util::sortCrsEntries;

  typedef typename Tpetra::CrsMatrix<Scalar,LO,GO,NT>::local_matrix_device_type
                   local_matrix_device_type;
  typedef typename local_matrix_device_type::StaticCrsGraphType graph_type;
  typedef typename graph_type::row_map_type::non_const_type rowptr_view_type;
  typedef typename graph_type::entries_type::non_const_type colind_view_type;
  typedef typename local_matrix_device_type::values_type::non_const_type 
                   vals_view_type;

  typedef typename rowptr_view_type::value_type index_type;
  typedef typename colind_view_type::value_type ordinal_type;
  typedef typename vals_view_type::value_type scalar_type;

  typedef typename std::vector<index_type> rowptr_type;
  typedef typename std::vector<ordinal_type> colind_type;
  typedef typename std::vector<scalar_type> vals_type;
  typedef typename colind_type::size_type size_type;

  // Map is not actually used, but is needed to instantiate Kokkos
  auto comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t INVALID =
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  auto dummy_map = Tpetra::Map<LO,GO,NT>(INVALID, 1, 0, comm);

  int max_num_entries_per_row = 7;  // should be odd
  int num_cols = 15; // should be odd

  rowptr_type  rowptr, rowptr2;
  colind_type  colind, colind2;
  vals_type vals, vals2;

  generate_crs_entries<scalar_type,ordinal_type,index_type>(
      rowptr, rowptr2, colind, colind2, vals, vals2, 
      max_num_entries_per_row, num_cols);

  {
    //
    // Sort the CRS entries by column index
    //

    size_type num_entries = colind.size();

    // Randomly shuffle values and column indices
    vals_type  vals_rand(num_entries);
    colind_type  colind_rand(num_entries);
    shuffle_crs_entries<scalar_type,ordinal_type,index_type>(
        colind_rand, vals_rand, rowptr, colind, vals);

    // Create mirror views of the CRS entries
    auto rowptr_views = 
         create_device_views<rowptr_view_type>(rowptr);
    auto colind_rand_views = 
         create_device_views<colind_view_type>(colind_rand);
    auto colind_rand_copy_views = 
         create_device_views<colind_view_type>(colind_rand);
    auto vals_rand_views = 
         create_device_views<vals_view_type>(vals_rand);

    //
    // Sort the GIDs and associated values
    //
    sortCrsEntries<rowptr_view_type,colind_view_type,vals_view_type>(
        rowptr_views.d, colind_rand_views.d, vals_rand_views.d);

    // Copy back to host
    rowptr_views.copy_to_host();
    colind_rand_views.copy_to_host();
    vals_rand_views.copy_to_host();

    // At this point, colind_rand and vals_rand should be sorted and equal to
    // their original arrays, respectively.
    TEST_COMPARE_ARRAYS(colind, colind_rand_views.h);

    // mfh 20 Mar 2018: This doesn't work for complex numbers.
    //
    //TEST_COMPARE_FLOATING_ARRAYS(vals, vals_rand_views.h, 1.e-12);
    {
      const bool equal_sizes =
        static_cast<std::size_t> (vals.size ()) == static_cast<std::size_t> (vals_rand_views.h.size ());
      TEST_ASSERT( equal_sizes );
      if (! equal_sizes) {
        for (std::size_t k = 0; k < vals.size (); ++k) {
          const auto diff = vals[k] - vals_rand_views.h[k];

          // FIXME (mfh 20 Mar 2018) Use a more appropriate tolerance for scalar_type.
          TEST_ASSERT( Kokkos::ArithTraits<scalar_type>::abs (diff) <= 1.e-12 );
        }
      }
    }

    //
    // Sort the GIDs w/o values
    //

    // Do the sort
    sortCrsEntries<rowptr_view_type, colind_view_type>(
        rowptr_views.d, colind_rand_copy_views.d);

    // Copy back to host
    rowptr_views.copy_to_host();
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
