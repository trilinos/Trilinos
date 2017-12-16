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

#include "KokkosGraph_GraphColor.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils.hpp"
#include "KokkosKernels_Handle.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using namespace KokkosGraph;
using namespace KokkosGraph::Experimental;

namespace Test {
template <typename crsMat_t, typename device>
int run_graphcolor(
    crsMat_t input_mat,
    ColoringAlgorithm coloring_algorithm,
    size_t &num_colors,
    typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type & vertex_colors){
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;


  typedef KokkosKernelsHandle
      <size_type,lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  kh.create_graph_coloring_handle(coloring_algorithm);


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();

  graph_color
    <KernelHandle,lno_view_t,lno_nnz_view_t> (&kh,num_rows_1, num_cols_1,
        input_mat.graph.row_map, input_mat.graph.entries);

  num_colors = kh.get_graph_coloring_handle()->get_num_colors();
  vertex_colors = kh.get_graph_coloring_handle()->get_vertex_colors();
  kh.destroy_graph_coloring_handle();
  return 0;
}

}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_coloring(lno_t numRows,size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::non_const_type   color_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  //typedef typename lno_view_t::non_const_value_type size_type;

  lno_t numCols = numRows;
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);

  typename lno_view_t::non_const_type sym_xadj;
  typename lno_nnz_view_t::non_const_type sym_adj;

  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<lno_view_t, lno_nnz_view_t,  typename lno_view_t::non_const_type, typename lno_nnz_view_t::non_const_type, device>
    (numRows, input_mat.graph.row_map, input_mat.graph.entries, sym_xadj, sym_adj);
  size_type numentries = sym_adj.dimension_0();
  scalar_view_t newValues("vals", numentries);

  graph_t static_graph (sym_adj, sym_xadj);
  input_mat = crsMat_t("CrsMatrix", numCols, newValues, static_graph);

  ColoringAlgorithm coloring_algorithms[] = {COLORING_DEFAULT, COLORING_SERIAL, COLORING_VB, COLORING_VBBIT, COLORING_VBCS, COLORING_EB};

  for (int ii = 0; ii < 6; ++ii){
    ColoringAlgorithm coloring_algorithm = coloring_algorithms[ii];
    color_view_t vector_colors;
    size_t num_colors;


    Kokkos::Impl::Timer timer1;
    crsMat_t output_mat;
    int res = run_graphcolor<crsMat_t, device>(input_mat, coloring_algorithm, num_colors, vector_colors);
    //double coloring_time = timer1.seconds();
    EXPECT_TRUE( (res == 0));


    const lno_t num_rows_1 = input_mat.numRows();
    const lno_t num_cols_1 = input_mat.numCols();
    lno_t num_conflict = KokkosKernels::Impl::kk_is_d1_coloring_valid
        <lno_view_t,lno_nnz_view_t, color_view_t, typename device::execution_space>
    (num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, vector_colors);

    lno_t conf = 0;
    {
      //also check the correctness of the validation code :)
      typename lno_view_t::HostMirror hrm = Kokkos::create_mirror_view (input_mat.graph.row_map);
      typename lno_nnz_view_t::HostMirror hentries = Kokkos::create_mirror_view (input_mat.graph.entries);
      typename color_view_t::HostMirror hcolor = Kokkos::create_mirror_view (vector_colors);
      Kokkos::deep_copy (hrm , input_mat.graph.row_map);
      Kokkos::deep_copy (hentries , input_mat.graph.entries);
      Kokkos::deep_copy (hcolor , vector_colors);

      for (lno_t i = 0; i < num_rows_1; ++i){
        const size_type b = hrm(i);
        const size_type e = hrm(i + 1);
        for (size_type j = b; j < e; ++j){
          lno_t d = hentries(j);
          if (i != d){
            if (hcolor(d) == hcolor(i)){
              conf++;
            }
          }
        }
      }
    }
    EXPECT_TRUE( (num_conflict == conf));

    EXPECT_TRUE( (num_conflict == 0));
  }
  //device::execution_space::finalize();

}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, graph ## _ ## graph_color ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_coloring<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 200, 10); \
  test_coloring<SCALAR,ORDINAL,OFFSET,DEVICE>(50000, 50000 * 30, 100, 10); \
}

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif
