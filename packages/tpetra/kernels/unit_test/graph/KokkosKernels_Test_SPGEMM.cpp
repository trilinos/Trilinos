
/*
//@HEADER
// ************************************************************************
//
//          KokkosKernels: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#  include "TpetraKernels_ETIHelperMacros.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <KokkosKernels_SPGEMM.hpp>
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils.hpp"
#include <Kokkos_Sparse_CrsMatrix.hpp>
#include <Kokkos_Concepts.hpp>
#include <string>
#include <stdexcept>

//const char *input_filename = "sherman1.mtx";
//const char *input_filename = "Si2.mtx";
//const char *input_filename = "wathen_30_30.mtx";
//const size_t expected_num_cols = 9906;



extern char * input_filename;

template <typename crsMat_t, typename device>
int run_spgemm(crsMat_t input_mat, crsMat_t input_mat2, KokkosKernels::Experimental::Graph::SPGEMMAlgorithm spgemm_algorithm, crsMat_t &result) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

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


  KokkosKernels::Experimental::Graph::spgemm_symbolic (
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
  KokkosKernels::Experimental::Graph::spgemm_numeric(
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


  KokkosKernels::Experimental::Util::kk_sort_graph<typename graph_t::row_map_type,
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

  if (nrows1 != nrows2) return false;
  if (nentries1 != nentries2) return false;
  if (nvals1 != nvals2) return false;

  KokkosKernels::Experimental::Util::kk_sort_graph
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
  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <typename graph_t::row_map_type, typename graph_t::row_map_type, typename lno_view_t::value_type,
      typename device::execution_space>(output_mat1.graph.row_map, output_mat2.graph.row_map, 0);
  if (!is_identical) return false;

  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(h_ent1, h_ent2, 0 );
  if (!is_identical) return false;

  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <scalar_view_t, scalar_view_t, typename scalar_view_t::value_type,
      typename device::execution_space>(h_vals1, h_vals2, 0.000001);
  if (!is_identical) return false;
  return true;
}

template <typename scalar_t, typename lno_t, typename device>
void test_spgemm(KokkosKernels::Experimental::Graph::SPGEMMAlgorithm spgemm_algorithm) {
  ASSERT_TRUE( (input_filename != NULL));
  device::execution_space::initialize();
  //device::execution_space::print_configuration(std::cout);

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;



  crsMat_t input_mat = KokkosKernels::Experimental::Util::read_kokkos_crst_matrix<crsMat_t>(input_filename);


  const int max_integer = 2147483647;
  std::string algo = "UNKNOWN";
  bool is_expected_to_fail = false;

  switch (spgemm_algorithm){
  case KokkosKernels::Experimental::Graph::SPGEMM_CUSPARSE:
    //TODO: add these test failure cases for cusparse too.
    algo = "SPGEMM_CUSPARSE";
#ifndef KERNELS_HAVE_CUSPARSE
    is_expected_to_fail = true;
#endif
    break;

  case KokkosKernels::Experimental::Graph::SPGEMM_MKL:
    algo = "SPGEMM_MKL";
#ifndef HAVE_TPETRAKERNELS_MKL
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
    if (input_mat.values.dimension_0() > max_integer){
      is_expected_to_fail = true;
    }

    if (!(Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device::memory_space>::accessible)){
      is_expected_to_fail = true;
    }
    break;

  case KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMSPEED:
    algo = "SPGEMM_KK_MEMSPEED";
    break;
  case KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED:
    algo = "SPGEMM_KK_SPEED";
    break;
  case KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY:
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

  double spgemm_time = timer1.seconds();
  if (!is_expected_to_fail){

    EXPECT_TRUE( (res == 0)) << algo;

    crsMat_t output_mat2;
    res = run_spgemm<crsMat_t, device>(input_mat, input_mat, KokkosKernels::Experimental::Graph::SPGEMM_DEBUG, output_mat2);


    bool is_identical = is_same_matrix<crsMat_t, device>(output_mat, output_mat2);

    EXPECT_TRUE(is_identical) << algo;

    //EXPECT_TRUE( equal) << algo;
  }
  device::execution_space::finalize();
}

#define INSTMACRO( SCALAR, LO, DEVICE) \
    test_spgemm<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY); \
    test_spgemm<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED); \
    test_spgemm<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMSPEED); \
    test_spgemm<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::SPGEMM_CUSPARSE); \
    test_spgemm<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::SPGEMM_MKL);




TEST (SPGEMM_TEST, SPGEMM) {
  TPETRAKERNELS_ETI_MANGLING_TYPEDEFS()
  TPETRAKERNELS_INSTANTIATE_SLD(INSTMACRO)
}
 

