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
#include <KokkosKernels_Handle.hpp>

#include <cstdlib>
#include <iostream>

#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <vector>

#include "KokkosGraph_graph_color.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_MyCRSMatrix.hpp"
#include "KokkosKernels_TestParameters.hpp"



void print_options(){

}
int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv){
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--threads" ) ) {
      params.use_threads = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--serial" ) ) {
      params.use_serial = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--openmp" ) ) {
      params.use_openmp = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--cuda" ) ) {
      params.use_cuda = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--repeat" ) ) {
      params.repeat = atoi( argv[++i] );
    }

    else if ( 0 == strcasecmp( argv[i] , "--chunksize" ) ) {
      params.chunk_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--teamsize" ) ) {
      params.team_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--vectorsize" ) ) {
      params.vector_size  = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--amtx" ) ) {
      params.a_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "--dynamic" ) ) {
      params.use_dynamic_scheduling = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--verbose" ) ) {
      params.verbose = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "--algorithm" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "COLORING_DEFAULT" ) ) {
        params.algorithm = 1;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_SERIAL" ) ) {
        params.algorithm = 2;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_VB" ) ) {
        params.algorithm = 3;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_VBBIT" ) ) {
        params.algorithm = 4;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_VBCS" ) ) {
        params.algorithm = 5;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_EB" ) ) {
        params.algorithm = 6;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_VBD" ) ) {
        params.algorithm = 7;
      }
      else if ( 0 == strcasecmp( argv[i] , "COLORING_VBDBIT" ) ) {
        params.algorithm = 8;
      }
      else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        print_options();
        return 1;
      }
    }
    else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      print_options();
      return 1;
    }
  }
  return 0;
}

namespace KokkosKernels{

namespace Experiment{


template <typename ExecSpace, typename crsGraph_t, typename crsGraph_t2 , typename crsGraph_t3 , typename TempMemSpace , typename PersistentMemSpace >
void run_experiment(
    crsGraph_t crsGraph, Parameters params){
  //using namespace KokkosSparse;
  using namespace KokkosGraph;
  using namespace KokkosGraph::Experimental;
  //using namespace KokkosSparse::Experimental;

  int algorithm = params.algorithm;
  int repeat = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize = params.shmemsize;
  int team_size = params.team_size;
  int use_dynamic_scheduling = params.use_dynamic_scheduling;
  int verbose = params.verbose;

  //char spgemm_step = params.spgemm_step;
  int vector_size = params.vector_size;

  typedef typename crsGraph_t3::row_map_type::non_const_type lno_view_t;
  typedef typename crsGraph_t3::entries_type::non_const_type lno_nnz_view_t;



  typedef typename lno_view_t::non_const_value_type size_type;
  typedef typename lno_nnz_view_t::non_const_value_type lno_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <size_type,lno_t, lno_t,
      ExecSpace, TempMemSpace,PersistentMemSpace > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);


  if (use_dynamic_scheduling){
    kh.set_dynamic_scheduling(true);
  }
  if (verbose){
    kh.set_verbose(true);
  }

  std::cout << "algorithm: " << algorithm << std::endl;

  for (int i = 0; i < repeat; ++i){

    switch (algorithm){
    case 1:
      kh.create_graph_coloring_handle(COLORING_DEFAULT);

      break;
    case 2:
      kh.create_graph_coloring_handle(COLORING_SERIAL);

      break;
    case 3:
      kh.create_graph_coloring_handle(COLORING_VB);
      break;
    case 4:
      kh.create_graph_coloring_handle(COLORING_VBBIT);

      break;
    case 5:
      kh.create_graph_coloring_handle(COLORING_VBCS);

      break;
    case 6:
      kh.create_graph_coloring_handle(COLORING_EB);
      break;

    case 7:
      kh.create_graph_coloring_handle(COLORING_VBD);
      break;

    case 8:
      kh.create_graph_coloring_handle(COLORING_VBDBIT);
      break;

    default:
      kh.create_graph_coloring_handle(COLORING_DEFAULT);

    }

    graph_color_symbolic(&kh,crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries);

    std::cout << std::endl <<
        "Time:" << kh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
        "Num colors:" << kh.get_graph_coloring_handle()->get_num_colors() << " "
        "Num Phases:" << kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
    std::cout << "\t"; KokkosKernels::Impl::print_1Dview(kh.get_graph_coloring_handle()->get_vertex_colors());

  }
}

template <typename size_type, typename lno_t,
          typename exec_space, typename hbm_mem_space, typename sbm_mem_space>
void run_multi_mem_experiment(Parameters params){

  typedef exec_space myExecSpace;
  typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
  typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

  typedef typename MyKokkosSparse::CrsMatrix<double, lno_t, myFastDevice, void, size_type > fast_crstmat_t;
  typedef typename fast_crstmat_t::StaticCrsGraphType fast_graph_t;
  //typedef typename fast_graph_t::row_map_type::non_const_type fast_row_map_view_t;
  //typedef typename fast_graph_t::entries_type::non_const_type   fast_cols_view_t;

  //typedef typename fast_graph_t::row_map_type::const_type const_fast_row_map_view_t;
  //typedef typename fast_graph_t::entries_type::const_type   const_fast_cols_view_t;

  typedef typename MyKokkosSparse::CrsMatrix<double, lno_t, mySlowExecSpace, void, size_type > slow_crstmat_t;
  typedef typename slow_crstmat_t::StaticCrsGraphType slow_graph_t;

  //typedef typename slow_graph_t::row_map_type::non_const_type slow_row_map_view_t;
  //typedef typename slow_graph_t::entries_type::non_const_type   slow_cols_view_t;
  //typedef typename slow_graph_t::row_map_type::const_type const_slow_row_map_view_t;
  //typedef typename slow_graph_t::entries_type::const_type   const_slow_cols_view_t;

  char *a_mat_file = params.a_mtx_bin_file;
  //char *b_mat_file = params.b_mtx_bin_file;
  //char *c_mat_file = params.c_mtx_bin_file;

  slow_graph_t a_slow_crsgraph, /*b_slow_crsgraph,*/ c_slow_crsgraph;
  fast_graph_t a_fast_crsgraph, /*b_fast_crsgraph,*/ c_fast_crsgraph;



  //read a and b matrices and store them on slow or fast memory.
  if (params.a_mem_space == 1){
    fast_crstmat_t a_fast_crsmat;
    a_fast_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
    a_fast_crsgraph = a_fast_crsmat.graph;
    a_fast_crsgraph.num_cols = a_fast_crsmat.numCols();

  }
  else {
    slow_crstmat_t a_slow_crsmat;
    a_slow_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
    a_slow_crsgraph = a_slow_crsmat.graph;
    a_slow_crsgraph.num_cols = a_slow_crsmat.numCols();
  }


  if (params.a_mem_space == 1){
    if (params.b_mem_space == 1){
      if (params.c_mem_space == 1){
        if (params.work_mem_space == 1){
           /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,fast_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else {
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,fast_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }

      }
      else {
        //C is in slow memory.
        if (params.work_mem_space == 1){
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,fast_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else {
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,fast_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
      }
    }
    else {
      //B is in slow memory
      if (params.c_mem_space == 1){
        if (params.work_mem_space == 1){
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,slow_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
                (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else {
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,slow_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
                (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }

      }
      else {
        //C is in slow memory.
        if (params.work_mem_space == 1){
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,slow_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
                (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else {
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,slow_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
                (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
      }

    }
  }
  else {
    //A is in slow memory
    if (params.b_mem_space == 1){
      if (params.c_mem_space == 1){
        if (params.work_mem_space == 1){
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,fast_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
                (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else {
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,fast_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
                (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }

      }
      else {
        //C is in slow memory.
        if (params.work_mem_space == 1){
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,fast_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
                (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else {
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,fast_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
                (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
      }
    }
    else {
      //B is in slow memory
      if (params.c_mem_space == 1){
        if (params.work_mem_space == 1){
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,slow_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
                (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else {
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,slow_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
                (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }

      }
      else {
        //C is in slow memory.
        if (params.work_mem_space == 1){
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,slow_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
                (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else {
          /*c_slow_crsgraph =*/
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, slow_graph_t,slow_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
                (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
      }

    }

  }
}



}
}

int main (int argc, char ** argv){

  typedef unsigned size_type;
  typedef int idx;
  //typedef int size_type;
  //typedef int idx;

  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs (params, argc, argv) ){
    return 1;
  }
  if (params.a_mtx_bin_file == NULL){
    std::cerr << "Provide a matrix file" << std::endl ;
    return 0;
  }
  std::cout << "Sizeof(idx):" << sizeof(idx) << " sizeof(size_type):" << sizeof(size_type) << std::endl;

  const int num_threads = params.use_openmp; // Assumption is that use_openmp variable is provided as number of threads
  const int device_id = 0;
  Kokkos::initialize( Kokkos::InitArguments( num_threads, -1, device_id ) );
  Kokkos::print_configuration(std::cout);

#if defined( KOKKOS_ENABLE_OPENMP )

  if (params.use_openmp) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::HostSpace>(
        params
        );
#else

    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>(
        params
        );
#endif
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  if (params.use_cuda) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::CudaHostPinnedSpace>(
        params
        );
#else
    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>(
        params
        );

#endif
  }

#endif

#if defined( KOKKOS_ENABLE_SERIAL )
  if (params.use_serial) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::Serial, Kokkos::Serial::memory_space, Kokkos::HostSpace>(
        params
        );
#else

    KokkosKernels::Experiment::run_multi_mem_experiment
    <size_type, idx, Kokkos::Serial, Kokkos::Serial::memory_space, Kokkos::Serial::memory_space>(
        params
        );
#endif
  }
#endif

  Kokkos::finalize();

  return 0;
}
