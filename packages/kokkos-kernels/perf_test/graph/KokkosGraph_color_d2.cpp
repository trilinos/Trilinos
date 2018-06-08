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
#include <string>


#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <vector>
#include "KokkosKernels_config.h"

#include "KokkosGraph_graph_color_d2.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_MyCRSMatrix.hpp"
#include "KokkosKernels_TestParameters.hpp"

using namespace KokkosGraph;
#ifdef KOKKOSKERNELS_INST_DOUBLE
typedef double kk_scalar_t;
#else
#ifdef KOKKOSKERNELS_INST_FLOAT
typedef float kk_scalar_t;
#endif
#endif

#ifdef KOKKOSKERNELS_INST_OFFSET_INT
typedef int kk_size_type;
#else
#ifdef KOKKOSKERNELS_INST_OFFSET_SIZE_T
typedef size_t kk_size_type;
#endif
#endif

#ifdef KOKKOSKERNELS_INST_ORDINAL_INT
typedef int kk_lno_t;
#else
#ifdef KOKKOSKERNELS_INST_ORDINAL_INT64_T
typedef int64_t kk_lno_t;
#endif
#endif


void print_options(std::ostream& os, const char* app_name, unsigned int indent=0)
{
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl 
     << spaces << app_name << " [parameters]" << std::endl 
     << std::endl
     << spaces << "Parameters:"
     << spaces << "  Parallelism (select one of the following):" << std::endl
     << spaces << "      serial <N>        Execute serially." << std::endl
     << spaces << "      threads <N>       Use N posix threads." << std::endl
     << spaces << "      openmp <N>        Use OpenMP with N threads." << std::endl
     << spaces << "      cuda              Use CUDA" << std::endl
     << std::endl
     << spaces << "  Required Parameters:" << std::endl 
     << spaces << "      amtx <filename>   Input file in Matrix Market format (.mtx)." << std::endl
     << std::endl
     << spaces << "      algorithm <algorithm_name>   Set the algorithm to use.  Allowable values are:" << std::endl
     << spaces << "                 COLORING_D2_MATRIX_SQUARED  - Distance-2 coloring using matrix-squared + Distance-1 coloring method." << std::endl
     << spaces << "                 COLORING_D2                 - Distance-2 coloring using traversal based method." << std::endl
     << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces << "      chunksize <N>     Set the chunk size." << std::endl
     << spaces << "      dynamic           Use dynamic scheduling." << std::endl
     << spaces << "      repeat <N>        Set number of test repetitions (Default: 6) " << std::endl
     << spaces << "      teamsize  <N>     Set the team size." << std::endl
     << spaces << "      vectorsize <N>    Set the vector size." << std::endl
     << spaces << "      verbose           Enable verbose mode (print timing + extra information" << std::endl
     << spaces << "      help              Print out command line help." << std::endl
     << spaces << " " << std::endl;
}


int parse_inputs (KokkosKernels::Experiment::Parameters &params, int argc, char **argv)
{
  bool got_required_param_amtx=false;
  bool got_required_param_algorithm=false;

  for(int i = 1; i < argc; ++i) 
  {
    if ( 0 == strcasecmp( argv[i] , "threads" ) ) 
    {
      params.use_threads = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "serial" ) ) 
    {
      params.use_serial = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) 
    {
      params.use_openmp = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) 
    {
      params.use_cuda = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "repeat" ) ) 
    {
      params.repeat = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "chunksize" ) ) 
    {
      params.chunk_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "teamsize" ) ) 
    {
      params.team_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "vectorsize" ) ) 
    {
      params.vector_size  = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "amtx" ) ) 
    {
      got_required_param_amtx = true;
      params.a_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "dynamic" ) ) 
    {
      params.use_dynamic_scheduling = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "verbose" ) ) 
    {
      params.verbose = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "algorithm" ) ) 
    {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "COLORING_D2_MATRIX_SQUARED" ) ) 
      {
        params.algorithm = 1;
        got_required_param_algorithm = true;
      }
      else if ( 0 == strcasecmp( argv[i], "COLORING_D2" ) )
      {
        params.algorithm = 2;
        got_required_param_algorithm = true;
      }
      else 
      {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        print_options(std::cout, argv[0]);
        return 1;
      }
    }
    else if ( 0 == strcasecmp( argv[i], "help") || 0 == strcasecmp(argv[i], "-h") )
    {
      print_options(std::cout, argv[0]);
      return 1;
    }
    else 
    {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      print_options(std::cout, argv[0]);
      return 1;
    }
  }

  if(!got_required_param_amtx)
  {
    std::cout << "Missing required parameter amtx" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if(!got_required_param_algorithm)
  {
    std::cout << "Missing required parameter algorithm" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if(!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda)
  {
    print_options(std::cout, argv[0]);
    return 1;
  }
  return 0;
}

namespace KokkosKernels {
namespace Experiment {

template <typename ExecSpace, typename crsGraph_t, typename crsGraph_t2 , typename crsGraph_t3 , typename TempMemSpace , typename PersistentMemSpace >
void run_experiment(crsGraph_t crsGraph, Parameters params)
{
  using namespace KokkosGraph;
  using namespace KokkosGraph::Experimental;

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

  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, kk_scalar_t, ExecSpace, TempMemSpace, PersistentMemSpace> KernelHandle;
    
  // Note: crsGraph.numRows() == number of vertices in the 'graph'
  //       crsGraph.entries.dimension_0() == number of edges in the 'graph'

  std::cout << "Num verts: " << crsGraph.numRows() << std::endl
            << "Num edges: " << crsGraph.entries.dimension_0() << std::endl;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);

  if (use_dynamic_scheduling)
  {
    kh.set_dynamic_scheduling(true);
  }

  if (verbose)
  {
    kh.set_verbose(true);
  }

  for (int i = 0; i < repeat; ++i)
  {

    switch (algorithm)
    {
    case 1:
      // kh.create_graph_coloring_handle(COLORING_SPGEMM);
      kh.create_graph_coloring_handle(COLORING_D2_MATRIX_SQUARED);
      break;
    case 2:
      kh.create_graph_coloring_handle(COLORING_D2);
      break;
    default:
      kh.create_graph_coloring_handle(COLORING_D2_MATRIX_SQUARED);
      break;
    }

    graph_color_d2(&kh, crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries);

    std::cout << "Time      : " << kh.get_graph_coloring_handle()->get_overall_coloring_time() << std::endl
              << "Num colors: " << kh.get_graph_coloring_handle()->get_num_colors() << std::endl
              << "Num Phases: " << kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
    std::cout << "\t"; 
    KokkosKernels::Impl::print_1Dview(kh.get_graph_coloring_handle()->get_vertex_colors());
    std::cout << std::endl;

  }
}



template <typename size_type, typename lno_t, typename exec_space, typename hbm_mem_space, typename sbm_mem_space>
void run_multi_mem_experiment(Parameters params)
{

  typedef exec_space myExecSpace;
  typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
  typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

  typedef typename MyKokkosSparse::CrsMatrix<double, lno_t, myFastDevice, void, size_type > fast_crstmat_t;
  typedef typename fast_crstmat_t::StaticCrsGraphType fast_graph_t;
  //typedef typename fast_graph_t::row_map_type::non_const_type fast_row_map_view_t;
  //typedef typename fast_graph_t::entries_type::non_const_type fast_cols_view_t;

  //typedef typename fast_graph_t::row_map_type::const_type const_fast_row_map_view_t;
  //typedef typename fast_graph_t::entries_type::const_type const_fast_cols_view_t;

  typedef typename MyKokkosSparse::CrsMatrix<double, lno_t, mySlowExecSpace, void, size_type > slow_crstmat_t;
  typedef typename slow_crstmat_t::StaticCrsGraphType slow_graph_t;

  //typedef typename slow_graph_t::row_map_type::non_const_type slow_row_map_view_t;
  //typedef typename slow_graph_t::entries_type::non_const_type slow_cols_view_t;
  //typedef typename slow_graph_t::row_map_type::const_type     const_slow_row_map_view_t;
  //typedef typename slow_graph_t::entries_type::const_type     const_slow_cols_view_t;

  char *a_mat_file = params.a_mtx_bin_file;

  slow_graph_t a_slow_crsgraph, /*b_slow_crsgraph,*/ c_slow_crsgraph;
  fast_graph_t a_fast_crsgraph, /*b_fast_crsgraph,*/ c_fast_crsgraph;


  //read a and b matrices and store them on slow or fast memory.
  if (params.a_mem_space == 1)
  {
    fast_crstmat_t a_fast_crsmat;
    a_fast_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
    a_fast_crsgraph = a_fast_crsmat.graph;
    a_fast_crsgraph.num_cols = a_fast_crsmat.numCols();

  }
  else 
  {
    slow_crstmat_t a_slow_crsmat;
    a_slow_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
    a_slow_crsgraph = a_slow_crsmat.graph;
    a_slow_crsgraph.num_cols = a_slow_crsmat.numCols();
  }


  if (params.a_mem_space == 1)
  {
    if (params.b_mem_space == 1)
    {
      if (params.c_mem_space == 1)
      {
        if (params.work_mem_space == 1)
        {
           /* c_fast_crsgraph = */
           KokkosKernels::Experiment::run_experiment
             <myExecSpace, fast_graph_t,fast_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
             (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else 
        {
          /* c_fast_crsgraph = */
              KokkosKernels::Experiment::run_experiment
                <myExecSpace, fast_graph_t,fast_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }

      }
      else 
      {
        //C is in slow memory.
        if (params.work_mem_space == 1)
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,fast_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
            (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else 
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,fast_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
            (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
        }
      }
    }
    else 
    {
      //B is in slow memory
      if (params.c_mem_space == 1)
      {
        if (params.work_mem_space == 1)
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,slow_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
            (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else 
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,slow_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
            (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }

      }
      else 
      {
        //C is in slow memory.
        if (params.work_mem_space == 1)
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,slow_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
            (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else 
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, fast_graph_t,slow_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
            (a_fast_crsgraph, /*b_slow_crsgraph,*/ params);
        }
      }

    }
  }
  else 
  {
    //A is in slow memory
    if (params.b_mem_space == 1)
    {
      if (params.c_mem_space == 1)
      {
        if (params.work_mem_space == 1)
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,fast_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
            (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else 
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,fast_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
            (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
      }
      else 
      {
        //C is in slow memory.
        if (params.work_mem_space == 1)
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,fast_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
            (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
        else 
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,fast_graph_t,slow_graph_t, sbm_mem_space, sbm_mem_space>
            (a_slow_crsgraph, /*b_fast_crsgraph,*/ params);
        }
      }
    }
    else 
    {
      //B is in slow memory
      if (params.c_mem_space == 1)
      {
        if (params.work_mem_space == 1)
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,slow_graph_t,fast_graph_t, hbm_mem_space, hbm_mem_space>
            (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else 
        {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,slow_graph_t,fast_graph_t, sbm_mem_space, sbm_mem_space>
            (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
      }
      else 
      {
        //C is in slow memory.
        if (params.work_mem_space == 1) 
        {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment
            <myExecSpace, slow_graph_t,slow_graph_t,slow_graph_t, hbm_mem_space, hbm_mem_space>
            (a_slow_crsgraph, /*b_slow_crsgraph,*/ params);
        }
        else 
        {
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



int main (int argc, char ** argv)
{
  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs (params, argc, argv) )
  {
    return 1;
  }

  if (params.a_mtx_bin_file == NULL)
  {
    std::cerr << "Provide a matrix file" << std::endl ;
    return 0;
  }

  std::cout << "Sizeof(kk_lno_t) : " << sizeof(kk_lno_t) << std::endl 
            << "Sizeof(size_type): " << sizeof(kk_size_type) << std::endl;

#if defined( KOKKOS_HAVE_OPENMP )

  if (params.use_openmp) 
  {
    Kokkos::OpenMP::initialize( params.use_openmp );
    Kokkos::OpenMP::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
        <kk_size_type, kk_lno_t, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::HostSpace>(params);
#else
    KokkosKernels::Experiment::run_multi_mem_experiment
        <kk_size_type, kk_lno_t, Kokkos::OpenMP, Kokkos::OpenMP::memory_space, Kokkos::OpenMP::memory_space>(params);
#endif
    Kokkos::OpenMP::finalize();
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  if (params.use_cuda) 
  {
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( 0 ) );
    Kokkos::Cuda::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
    <kk_size_type, kk_lno_t, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::CudaHostPinnedSpace>(params);
#else
    KokkosKernels::Experiment::run_multi_mem_experiment
    <kk_size_type, kk_lno_t, Kokkos::Cuda, Kokkos::Cuda::memory_space, Kokkos::Cuda::memory_space>(params);
#endif
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_SERIAL )
  if (params.use_serial) 
  {
    Kokkos::Serial::initialize( params.use_openmp );
    Kokkos::Serial::print_configuration(std::cout);
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment
    <kk_size_type, kk_lno_t, Kokkos::Serial, Kokkos::Serial::memory_space, Kokkos::HostSpace>(params);
#else
    KokkosKernels::Experiment::run_multi_mem_experiment
    <kk_size_type, kk_lno_t, Kokkos::Serial, Kokkos::Serial::memory_space, Kokkos::Serial::memory_space>(params);
#endif
    Kokkos::Serial::finalize();
  }
#endif

  return 0;
}

