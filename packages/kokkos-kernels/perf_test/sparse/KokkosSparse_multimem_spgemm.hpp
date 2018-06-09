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

#include "KokkosKernels_MyCRSMatrix.hpp"
#include "KokkosSparse_run_spgemm.hpp"
namespace KokkosKernels{

namespace Experiment{

  template <typename size_type, typename lno_t, typename scalar_t,
            typename exec_space, typename hbm_mem_space, typename sbm_mem_space>
  void run_multi_mem_spgemm(Parameters params){

    typedef exec_space myExecSpace;
    typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
    typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

    typedef typename MyKokkosSparse::CrsMatrix<scalar_t, lno_t, myFastDevice, void, size_type > fast_crstmat_t;
    //typedef typename fast_crstmat_t::StaticCrsGraphType fast_graph_t;
    //typedef typename fast_crstmat_t::row_map_type::non_const_type fast_row_map_view_t;
    typedef typename fast_crstmat_t::index_type::non_const_type   fast_cols_view_t;
    typedef typename fast_crstmat_t::values_type::non_const_type fast_values_view_t;
    typedef typename fast_crstmat_t::row_map_type::const_type const_fast_row_map_view_t;
    typedef typename fast_crstmat_t::index_type::const_type   const_fast_cols_view_t;
    typedef typename fast_crstmat_t::values_type::const_type const_fast_values_view_t;

    typedef typename MyKokkosSparse::CrsMatrix<scalar_t, lno_t, mySlowExecSpace, void, size_type > slow_crstmat_t;
    //typedef typename slow_crstmat_t::StaticCrsGraphType slow_graph_t;
    //typedef typename slow_crstmat_t::row_map_type::non_const_type slow_row_map_view_t;
    typedef typename slow_crstmat_t::index_type::non_const_type   slow_cols_view_t;
    typedef typename slow_crstmat_t::values_type::non_const_type slow_values_view_t;
    typedef typename slow_crstmat_t::row_map_type::const_type const_slow_row_map_view_t;
    typedef typename slow_crstmat_t::index_type::const_type   const_slow_cols_view_t;
    typedef typename slow_crstmat_t::values_type::const_type const_slow_values_view_t;

    char *a_mat_file = params.a_mtx_bin_file;
    char *b_mat_file = params.b_mtx_bin_file;
    char *c_mat_file = params.c_mtx_bin_file;

    slow_crstmat_t a_slow_crsmat, b_slow_crsmat, c_slow_crsmat;
    fast_crstmat_t a_fast_crsmat, b_fast_crsmat, c_fast_crsmat;

    //read a and b matrices and store them on slow or fast memory.

    if (params.a_mem_space == 1){
      a_fast_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
    }
    else {
      a_slow_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
    }


    if ((b_mat_file == NULL || strcmp(b_mat_file, a_mat_file) == 0) && params.b_mem_space == params.a_mem_space){
      std::cout << "Using A matrix for B as well" << std::endl;
      b_fast_crsmat = a_fast_crsmat;
      b_slow_crsmat = a_slow_crsmat;
    }
    else if (params.b_mem_space == 1){
      if (b_mat_file == NULL) b_mat_file = a_mat_file;
      b_fast_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(b_mat_file);
    }
    else {
      if (b_mat_file == NULL) b_mat_file = a_mat_file;
      b_slow_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(b_mat_file);
    }

    if (params.a_mem_space == 1){
      if (params.b_mem_space == 1){
        if (params.c_mem_space == 1){
          if (params.work_mem_space == 1){
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,fast_crstmat_t,fast_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_fast_crsmat, b_fast_crsmat, params);
          }
          else {
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,fast_crstmat_t,fast_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_fast_crsmat, b_fast_crsmat, params);
          }

        }
        else {
          //C is in slow memory.
          if (params.work_mem_space == 1){
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,fast_crstmat_t,slow_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_fast_crsmat, b_fast_crsmat, params);
          }
          else {
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,fast_crstmat_t,slow_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_fast_crsmat, b_fast_crsmat, params);
          }
        }
      }
      else {
        //B is in slow memory
        if (params.c_mem_space == 1){
          if (params.work_mem_space == 1){
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,slow_crstmat_t,fast_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_fast_crsmat, b_slow_crsmat, params);
          }
          else {
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,slow_crstmat_t,fast_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_fast_crsmat, b_slow_crsmat, params);
          }

        }
        else {
          //C is in slow memory.
          if (params.work_mem_space == 1){
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,slow_crstmat_t,slow_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_fast_crsmat, b_slow_crsmat, params);
          }
          else {
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, fast_crstmat_t,slow_crstmat_t,slow_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_fast_crsmat, b_slow_crsmat, params);
          }
        }

      }
    }
    else {
      //A is in slow memory
      if (params.b_mem_space == 1){
        if (params.c_mem_space == 1){
          if (params.work_mem_space == 1){
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,fast_crstmat_t,fast_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_slow_crsmat, b_fast_crsmat, params);
          }
          else {
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,fast_crstmat_t,fast_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_slow_crsmat, b_fast_crsmat, params);
          }

        }
        else {
          //C is in slow memory.
          if (params.work_mem_space == 1){
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,fast_crstmat_t,slow_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_slow_crsmat, b_fast_crsmat, params);
          }
          else {
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,fast_crstmat_t,slow_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_slow_crsmat, b_fast_crsmat, params);
          }
        }
      }
      else {
        //B is in slow memory
        if (params.c_mem_space == 1){
          if (params.work_mem_space == 1){
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,slow_crstmat_t,fast_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_slow_crsmat, b_slow_crsmat, params);
          }
          else {
            c_fast_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,slow_crstmat_t,fast_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_slow_crsmat, b_slow_crsmat, params);
          }

        }
        else {
          //C is in slow memory.
          if (params.work_mem_space == 1){
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,slow_crstmat_t,slow_crstmat_t, hbm_mem_space, hbm_mem_space>
                  (a_slow_crsmat, b_slow_crsmat, params);
          }
          else {
            c_slow_crsmat =
                KokkosKernels::Experiment::run_experiment
                  <myExecSpace, slow_crstmat_t,slow_crstmat_t,slow_crstmat_t, sbm_mem_space, sbm_mem_space>
                  (a_slow_crsmat, b_slow_crsmat, params);
          }
        }

      }

    }


    if (c_mat_file != NULL){
      if (params.c_mem_space == 1){

        fast_cols_view_t sorted_adj("sorted adj", c_fast_crsmat.graph.entries.extent(0));
        fast_values_view_t sorted_vals("sorted vals", c_fast_crsmat.graph.entries.extent(0));

        KokkosKernels::Impl::kk_sort_graph
        <const_fast_row_map_view_t, const_fast_cols_view_t, const_fast_values_view_t, fast_cols_view_t, fast_values_view_t, myExecSpace>(
            c_fast_crsmat.graph.row_map,
            c_fast_crsmat.graph.entries,
            c_fast_crsmat.values, sorted_adj, sorted_vals);

        KokkosKernels::Impl::write_graph_bin(
            (lno_t) (c_fast_crsmat.numRows()),
            (size_type) (c_fast_crsmat.graph.entries.extent(0)),
            c_fast_crsmat.graph.row_map.data(),
            sorted_adj.data(),
            sorted_vals.data(),
            c_mat_file);
      }
      else {
        slow_cols_view_t sorted_adj("sorted adj", c_fast_crsmat.graph.entries.extent(0));
        slow_values_view_t sorted_vals("sorted vals", c_fast_crsmat.graph.entries.extent(0));

        KokkosKernels::Impl::kk_sort_graph<
        const_slow_row_map_view_t, const_slow_cols_view_t, const_slow_values_view_t, slow_cols_view_t, slow_values_view_t, myExecSpace>(
            c_slow_crsmat.graph.row_map,
            c_slow_crsmat.graph.entries,
            c_slow_crsmat.values, sorted_adj, sorted_vals);

        KokkosKernels::Impl::write_graph_bin(
            (lno_t) c_slow_crsmat.numRows(),
            (size_type) c_slow_crsmat.graph.entries.extent(0),
            c_slow_crsmat.graph.row_map.data(),
            sorted_adj.data(),
            sorted_vals.data(),
                    c_mat_file);
      }
    }
  }


}
}
