#include<gtest/gtest.h>   
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<Kokkos_ArithTraits.hpp>

#include<KokkosSparse_CrsMatrix.hpp>
#include<KokkosSparse_spadd.hpp>
#include<KokkosKernels_TestUtils.hpp>
#include<KokkosKernels_IOUtils.hpp>
#include<KokkosKernels_Utils.hpp>

#include<algorithm>   //for std::random_shuffle
#include<random>      // for std::default_random_engine
#include<cstdlib>     //for rand
#include<type_traits> //for std::is_same

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

//Create a random square matrix for testing mat-mat addition kernels
template <typename crsMat_t, typename ordinal_type>
crsMat_t randomMatrix(ordinal_type nrows, ordinal_type minNNZ, ordinal_type maxNNZ, bool sortRows)
{
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type size_type_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename size_type_view_t::non_const_value_type size_type;  //rowptr type
  typedef typename lno_view_t::non_const_value_type lno_t;            //colind type
  static_assert(std::is_same<ordinal_type, lno_t>::value, "ordinal_type should be same as lno_t from crsMat_t");
  //first, populate rowmap
  size_type_view_t rowmap("rowmap", nrows + 1);
  typename size_type_view_t::HostMirror h_rowmap = Kokkos::create_mirror_view(rowmap);
  size_type nnz = 0;
  for(lno_t i = 0; i < nrows; i++)
  {
    size_type rowEntries = rand() % (maxNNZ - minNNZ + 1) + minNNZ;
    h_rowmap(i) = nnz;
    nnz += rowEntries;
  }
  h_rowmap(nrows) = nnz;
  Kokkos::deep_copy(rowmap, h_rowmap);
  //allocate values and entries
  scalar_view_t values("values", nnz);
  //populate values
  typename scalar_view_t::HostMirror h_values = Kokkos::create_mirror_view(values);
  for(size_type i = 0; i < nnz; i++)
  {
    h_values(i) = rand() / RAND_MAX;
  }
  Kokkos::deep_copy(values, h_values);
  //populate entries (make sure no entry is repeated within a row)
  lno_view_t entries("entries", nnz);
  typename lno_view_t::HostMirror h_entries = Kokkos::create_mirror_view(entries);
  std::vector<lno_t> indices(nrows);
  auto re = std::default_random_engine(0);
  for(lno_t i = 0; i < nrows; i++)
  {
    for(lno_t j = 0; j < nrows; j++)
    {
      indices[j] = j;
    }
    std::shuffle(indices.begin(), indices.end(), re);
    size_type rowStart = h_rowmap(i);
    size_type rowCount = h_rowmap(i + 1) - rowStart;
    if(sortRows)
    {
      std::sort(indices.begin(), indices.begin() + rowCount);
    }
    for(size_type j = 0; j < rowCount; j++)
    {
      h_entries(rowStart + j) = indices[j];
    }
  }
  Kokkos::deep_copy(entries, h_entries);
  return crsMat_t("test matrix", nrows, nrows, nnz, values, rowmap, entries);
}

template <typename scalar_t, typename lno_t, typename size_type, class Device>
void test_spadd(lno_t numRows, size_type minNNZ, size_type maxNNZ, bool sortRows)
{
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device, void, size_type> crsMat_t;

  typedef typename crsMat_t::row_map_type::non_const_type row_map_type;
  typedef typename crsMat_t::index_type::non_const_type entries_type;
  typedef typename crsMat_t::values_type::non_const_type values_type;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t,
  typename Device::execution_space, typename Device::memory_space, typename Device::memory_space> KernelHandle;

  KernelHandle handle;
  handle.create_spadd_handle(sortRows);
  crsMat_t A = randomMatrix<crsMat_t, lno_t>(numRows, minNNZ, maxNNZ, sortRows);
  crsMat_t B = randomMatrix<crsMat_t, lno_t>(numRows, minNNZ, maxNNZ, sortRows);
  //Matrices from randomMatrix are always square
  lno_t numCols = numRows;
  row_map_type c_row_map("C row map", numRows + 1);
  auto addHandle = handle.get_spadd_handle();
  KokkosSparse::Experimental::spadd_symbolic<
    KernelHandle,
    typename row_map_type::const_type,
    typename entries_type::const_type,
    typename row_map_type::const_type,
    typename entries_type::const_type,
    row_map_type,
    entries_type>
  (&handle, A.graph.row_map, A.graph.entries, B.graph.row_map, B.graph.entries, c_row_map);
  values_type c_values("C values", addHandle->get_max_result_nnz());
  entries_type c_entries("C entries", addHandle->get_max_result_nnz());
  KokkosSparse::Experimental::spadd_numeric<
    KernelHandle,
    typename row_map_type::const_type,
    typename entries_type::const_type,
    scalar_t, typename values_type::const_type,
    typename row_map_type::const_type,
    typename entries_type::const_type,
    scalar_t, typename values_type::const_type,
    row_map_type, entries_type, values_type>
    (&handle, A.graph.row_map, A.graph.entries, A.values, 1,
     B.graph.row_map, B.graph.entries, B.values, 1,
     c_row_map, c_entries, c_values);
  //done with handle
  //create C using CRS arrays
  crsMat_t C("C", numRows, numRows, addHandle->get_max_result_nnz(), c_values, c_row_map, c_entries);
  handle.destroy_spadd_handle();
  auto Avalues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  auto Arowmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto Aentries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto Bvalues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B.values);
  auto Browmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B.graph.row_map);
  auto Bentries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B.graph.entries);
  auto Cvalues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C.values);
  auto Crowmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C.graph.row_map);
  auto Centries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C.graph.entries);
  using KAT = Kokkos::ArithTraits<scalar_t>;
  auto zero = KAT::zero();
  auto eps = KAT::epsilon();
  //check that C is correct and sorted, row-by-row
  for(lno_t row = 0; row < numRows; row++)
  {
    std::vector<scalar_t> correct(numCols, zero);
    std::vector<bool> nonzeros(numCols, false);
    for(size_type i = Arowmap(row); i < Arowmap(row + 1); i++)
    {
      correct[Aentries(i)] += Avalues(i);
      nonzeros[Aentries(i)] = true;
    }
    for(size_type i = Browmap(row); i < Browmap(row + 1); i++)
    {
      correct[Bentries(i)] += Bvalues(i);
      nonzeros[Bentries(i)] = true;
    }
    size_type nz = 0;
    for(lno_t i = 0; i < numCols; i++)
    {
      if(nonzeros[i])
        nz++;
    }
    //make sure C has the right number of entries
    auto actualNZ = Crowmap(row + 1) - Crowmap(row);
    ASSERT_EQ(actualNZ, nz) << "A+B row " << row << " has " << actualNZ << " entries but should have " << nz;
    //make sure C's indices are sorted
    for(size_type i = Crowmap(row) + 1; i < Crowmap(row + 1); i++)
    {
      ASSERT_LE(Centries(i - 1), Centries(i)) << "C row " << row << " is not sorted";
    }
    //make sure C has the correct values
    for(size_type i = Crowmap(row); i < Crowmap(row + 1); i++)
    {
      scalar_t Cval = Cvalues(i);
      lno_t Ccol = Centries(i);
      //Check that result is correct to 1 ULP
      ASSERT_LE(KAT::abs(correct[Ccol] - Cval), KAT::abs(correct[Ccol] * eps)) << "A+B row " << row << ", column " << Ccol << " has value " << Cval << " but should be " << correct[Ccol];
    }
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory,sparse ## _ ## spadd ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_spadd<SCALAR,ORDINAL,OFFSET,DEVICE> (10, 0, 3, true); \
  test_spadd<SCALAR,ORDINAL,OFFSET,DEVICE> (10, 0, 3, false); \
  test_spadd<SCALAR,ORDINAL,OFFSET,DEVICE> (100, 50, 100, true); \
  test_spadd<SCALAR,ORDINAL,OFFSET,DEVICE> (100, 50, 100, false); \
}

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

