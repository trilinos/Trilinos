/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosKernels_SparseUtils.hpp"
#include <Kokkos_Concepts.hpp>
#include <string>
#include <stdexcept>

#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>


#include<KokkosKernels_IOUtils.hpp>


//const char *input_filename = "sherman1.mtx";
//const char *input_filename = "Si2.mtx";
//const char *input_filename = "wathen_30_30.mtx";
//const size_t expected_num_cols = 9906;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test {

template <typename crsMat_t, typename device>
int run_spgemm(crsMat_t input_mat, crsMat_t input_mat2, KokkosSparse::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &result) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;



  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <size_type,lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  //kh.set_verbose(true);

  kh.create_spgemm_handle(spgemm_algorithm);


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_rows_2 = input_mat2.numRows();
  const size_t num_cols_2 = input_mat2.numCols();

  const size_t num_cols_1 = input_mat.numCols();
  bool equal = num_rows_2 == num_cols_1;
  if (!equal) return 1;



  lno_view_t row_mapC ("non_const_lnow_row", num_rows_1 + 1);
  lno_nnz_view_t  entriesC;
  scalar_view_t valuesC;


  spgemm_symbolic (
      &kh,
      num_rows_1,
      num_rows_2,
      num_cols_2,
      input_mat.graph.row_map,
      input_mat.graph.entries,
      false,
      input_mat2.graph.row_map,
      input_mat2.graph.entries,
      false,
      row_mapC
  );

  size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size){
    entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
    valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
  }
  spgemm_numeric(
      &kh,
      num_rows_1,
      num_rows_2,
      num_cols_2,
      input_mat.graph.row_map,
      input_mat.graph.entries,
      input_mat.values,
      false,

      input_mat2.graph.row_map,
      input_mat2.graph.entries,
      input_mat2.values,
      false,
      row_mapC,
      entriesC,
      valuesC
  );


  graph_t static_graph (entriesC, row_mapC);
  crsMat_t crsmat("CrsMatrix", num_cols_2, valuesC, static_graph);
  result = crsmat;
  kh.destroy_spgemm_handle();

  return 0;
}
template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat1, crsMat_t output_mat2){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows1 = output_mat1.graph.row_map.dimension_0();
  size_t nentries1 = output_mat1.graph.entries.dimension_0() ;
  size_t nvals1 = output_mat1.values.dimension_0();

  size_t nrows2 = output_mat2.graph.row_map.dimension_0();
  size_t nentries2 = output_mat2.graph.entries.dimension_0() ;
  size_t nvals2 = output_mat2.values.dimension_0();


  lno_nnz_view_t h_ent1 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries1);
  scalar_view_t h_vals1 (Kokkos::ViewAllocateWithoutInitializing("v1"), nvals1);


  KokkosKernels::Impl::kk_sort_graph<typename graph_t::row_map_type,
    typename graph_t::entries_type,
    typename crsMat_t::values_type,
    lno_nnz_view_t,
    scalar_view_t,
    typename device::execution_space
    >(
    output_mat1.graph.row_map, output_mat1.graph.entries, output_mat1.values,
    h_ent1, h_vals1
  );

  lno_nnz_view_t h_ent2 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries2);
  scalar_view_t h_vals2 (Kokkos::ViewAllocateWithoutInitializing("v1"), nvals2);

  if (nrows1 != nrows2) { 
     std::cout << "nrows1:" << nrows1 << " nrows2:" << nrows2 << std::endl;
     return false;
  }
  if (nentries1 != nentries2) {
    std::cout << "nentries1:" << nentries1 << " nentries2:" << nentries2 << std::endl;
    return false;
  }
  if (nvals1 != nvals2) {
    std::cout << "nvals1:" << nvals1 << " nvals2:" << nvals2 << std::endl;
    return false;
  }

  KokkosKernels::Impl::kk_sort_graph
      <typename graph_t::row_map_type,
      typename graph_t::entries_type,
      typename crsMat_t::values_type,
      lno_nnz_view_t,
      scalar_view_t,
      typename device::execution_space
      >(
      output_mat2.graph.row_map, output_mat2.graph.entries, output_mat2.values,
      h_ent2, h_vals2
    );

  bool is_identical = true;
  is_identical = KokkosKernels::Impl::kk_is_identical_view
      <typename graph_t::row_map_type, typename graph_t::row_map_type, typename lno_view_t::value_type,
      typename device::execution_space>(output_mat1.graph.row_map, output_mat2.graph.row_map, 0);
  //KokkosKernels::Impl::kk_print_1Dview(output_mat2.graph.row_map);

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat1.graph.row_map);
    KokkosKernels::Impl::kk_print_1Dview(output_mat2.graph.row_map);
    return false;
  }

  is_identical = KokkosKernels::Impl::kk_is_identical_view
      <lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(h_ent1, h_ent2, 0 );

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(h_ent1);
    KokkosKernels::Impl::kk_print_1Dview(h_ent2);
    return false;
  }


  typedef typename Kokkos::Details::ArithTraits<typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type,float>::value?2*1e-3:1e-7;


  is_identical = KokkosKernels::Impl::kk_is_identical_view
      <scalar_view_t, scalar_view_t, eps_type,
      typename device::execution_space>(h_vals1, h_vals2, eps);

  if (!is_identical) {
    std::cout << "values are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat1.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat2.values);

    return false;
  }
  return true;
}
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spgemm(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {

  using namespace Test;
  //device::execution_space::initialize();
  //device::execution_space::print_configuration(std::cout);

  typedef CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
  //typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  //typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  //typedef typename crsMat_t::values_type::non_const_type scalar_view_t;


  lno_t numCols = numRows;
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);


  crsMat_t output_mat2;
  run_spgemm<crsMat_t, device>(input_mat, input_mat, SPGEMM_DEBUG, output_mat2);

  SPGEMMAlgorithm algorithms [] = {SPGEMM_KK_MEMORY, SPGEMM_KK_SPEED, SPGEMM_KK_MEMSPEED, SPGEMM_CUSPARSE,SPGEMM_MKL};

  for (int ii = 0; ii < 5; ++ii){

    SPGEMMAlgorithm spgemm_algorithm = algorithms[ii];

    const uint64_t max_integer = 2147483647;
    std::string algo = "UNKNOWN";
    bool is_expected_to_fail = false;

    switch (spgemm_algorithm){
    case SPGEMM_CUSPARSE:
      //TODO: add these test failure cases for cusparse too.
      algo = "SPGEMM_CUSPARSE";
#ifndef KERNELS_HAVE_CUSPARSE
      is_expected_to_fail = true;
#endif
      break;

    case SPGEMM_MKL:
      algo = "SPGEMM_MKL";
#ifndef HAVE_KOKKOSKERNELS_MKL
      is_expected_to_fail = true;
#endif
      //MKL requires scalar to be either float or double
      if (!(Kokkos::Impl::is_same<float,scalar_t>::value || Kokkos::Impl::is_same<double,scalar_t>::value)){
        is_expected_to_fail = true;
      }
      //mkl requires local ordinals to be int.
      if (!(Kokkos::Impl::is_same<int,lno_t>::value)){
        is_expected_to_fail = true;
      }
      //if size_type is larger than int, mkl casts it to int.
      //it will fail if casting cause overflow.
      if (input_mat.values.extent(0) > max_integer){
        is_expected_to_fail = true;
      }

      if (!(Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device::memory_space>::accessible)){
        is_expected_to_fail = true;
      }
      break;

    case SPGEMM_KK_MEMSPEED:
      algo = "SPGEMM_KK_MEMSPEED";
      break;
    case SPGEMM_KK_SPEED:
      algo = "SPGEMM_KK_SPEED";
      break;
    case SPGEMM_KK_MEMORY:
      algo = "SPGEMM_KK_MEMORY";
      break;
    default:
      break;
    }

    Kokkos::Impl::Timer timer1;
    crsMat_t output_mat;

    bool failed = false;
    int res = 0;
    try{
      res = run_spgemm<crsMat_t, device>(input_mat, input_mat, spgemm_algorithm, output_mat);
    }
    catch (const char *message){
      EXPECT_TRUE(is_expected_to_fail) << algo;
      failed = true;
    }
    catch (std::string message){
      EXPECT_TRUE(is_expected_to_fail)<< algo;
      failed = true;
    }
    catch (std::exception& e){
      EXPECT_TRUE(is_expected_to_fail)<< algo;
      failed = true;
    }
    EXPECT_TRUE((failed == is_expected_to_fail));

    //double spgemm_time = timer1.seconds();

    timer1.reset();
    if (!is_expected_to_fail){

      EXPECT_TRUE( (res == 0)) << algo;
      bool is_identical = is_same_matrix<crsMat_t, device>(output_mat, output_mat2);
      EXPECT_TRUE(is_identical) << algo;
      //EXPECT_TRUE( equal) << algo;
    }
    //std::cout << "algo:" << algo << " spgemm_time:" << spgemm_time << " output_check_time:" << timer1.seconds() << std::endl;
  }
  //device::execution_space::finalize();
}



#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ## spgemm ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(10000, 10000 * 30, 500, 10); \
}

//test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 100, 10);
//test_spgemm<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 200, 10);

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


