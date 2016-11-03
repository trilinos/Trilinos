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

#include <KokkosKernels_GraphColor.hpp>
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils.hpp"
#include <Kokkos_Sparse_CrsMatrix.hpp>
#include "KokkosKernels_Handle.hpp"

//const char *input_filename = "sherman1.mtx";
//const char *input_filename = "Si2.mtx";
//const char *input_filename = "wathen_30_30.mtx";
extern char * input_filename;

template <typename crsMat_t, typename device>
int run_graphcolor(
    crsMat_t input_mat,
    KokkosKernels::Experimental::Graph::ColoringAlgorithm coloring_algorithm,
    size_t &num_colors,
    typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type & vertex_colors){
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  kh.create_graph_coloring_handle(coloring_algorithm);


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();

  KokkosKernels::Experimental::Graph::graph_color_symbolic
    <KernelHandle,lno_view_t,lno_nnz_view_t> (&kh,num_rows_1, num_cols_1,
        input_mat.graph.row_map, input_mat.graph.entries);

  num_colors = kh.get_graph_coloring_handle()->get_num_colors();
  vertex_colors = kh.get_graph_coloring_handle()->get_vertex_colors();
  kh.destroy_graph_coloring_handle();
  return 0;
}


template <typename scalar_t, typename lno_t, typename device>
void test_coloring(KokkosKernels::Experimental::Graph::ColoringAlgorithm coloring_algorithm) {
  ASSERT_TRUE( (input_filename != NULL));

  device::execution_space::initialize();
  //device::execution_space::print_configuration(std::cout);

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::non_const_type   color_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename lno_view_t::non_const_value_type size_type;

  crsMat_t input_mat = KokkosKernels::Experimental::Util::read_kokkos_crst_matrix<crsMat_t>(input_filename);


  color_view_t vector_colors;
  size_t num_colors;


  Kokkos::Impl::Timer timer1;
  crsMat_t output_mat;
  int res = run_graphcolor<crsMat_t, device>(input_mat, coloring_algorithm, num_colors, vector_colors);
  double coloring_time = timer1.seconds();
  EXPECT_TRUE( (res == 0));


  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();
  size_t num_conflict = KokkosKernels::Experimental::Util::kk_is_d1_coloring_valid
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

    for (size_t i = 0; i < num_rows_1; ++i){
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

  device::execution_space::finalize();

}

#define INSTMACRO( SCALAR, LO, DEVICE) \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_DEFAULT); \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_SERIAL); \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_VB); \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_VBBIT); \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_VBCS); \
    test_coloring<SCALAR,LO,DEVICE>(KokkosKernels::Experimental::Graph::COLORING_EB);



TEST (COLORING_TEST, GC) {
  TPETRAKERNELS_ETI_MANGLING_TYPEDEFS()
  TPETRAKERNELS_INSTANTIATE_SLD(INSTMACRO)
}
 

