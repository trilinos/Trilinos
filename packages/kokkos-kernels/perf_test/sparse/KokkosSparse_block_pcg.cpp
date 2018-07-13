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

#include <KokkosKernels_config.h>
#if defined(KOKKOSKERNELS_INST_DOUBLE) &&  \
    defined(KOKKOSKERNELS_INST_OFFSET_INT) && \
    defined(KOKKOSKERNELS_INST_ORDINAL_INT)
#include "KokkosSparse_pcg.hpp"

#include "KokkosKernels_Utils.hpp"
#include <iostream>
#include "KokkosKernels_IOUtils.hpp"

#define MAXVAL 1

#define SIZE_TYPE int
#define INDEX_TYPE int
#define SCALAR_TYPE double
unsigned cg_iteration_limit = 10;



template<typename scalar_view_t>
scalar_view_t create_x_vector(INDEX_TYPE nv, SCALAR_TYPE max_value = 1.0){
  scalar_view_t kok_x ("X", nv);

  typename scalar_view_t::HostMirror h_x =  Kokkos::create_mirror_view (kok_x);


  for (INDEX_TYPE i = 0; i < nv; ++i){
    SCALAR_TYPE r = static_cast <SCALAR_TYPE> (rand()) / static_cast <SCALAR_TYPE> (RAND_MAX / max_value);
    h_x(i) = r;
  }
  Kokkos::deep_copy (kok_x, h_x);
  return kok_x;
}



template <typename crsMat_t, typename vector_t>
vector_t create_y_vector(crsMat_t crsMat, vector_t x_vector){
  vector_t y_vector ("Y VECTOR", crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 1, y_vector);
  return y_vector;
}

template <typename ExecSpace, typename crsMat_t>
void run_point_experiment(
    crsMat_t crsmat,
	typename  crsMat_t::values_type::non_const_type kok_x_original){


	  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
	  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
	  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
	  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

	  typedef typename lno_nnz_view_t::value_type lno_t;
	  typedef typename lno_view_t::value_type size_type;
	  typedef typename scalar_view_t::value_type scalar_t;
	  INDEX_TYPE nv = crsmat.numRows();
	  //scalar_view_t kok_x_original = create_x_vector<scalar_view_t>(nv, MAXVAL);

	  //KokkosKernels::Impl::print_1Dview(kok_x_original);
	  scalar_view_t kok_b_vector = create_y_vector(crsmat, kok_x_original);

	  //create X vector
	  scalar_view_t kok_x_vector("kok_x_vector", nv);


	  double solve_time = 0;
	  const double   cg_iteration_tolerance     = 1e-7 ;

	  KokkosKernels::Experimental::Example::CGSolveResult cg_result ;




	  typedef KokkosKernels::Experimental::KokkosKernelsHandle
	        < size_type,
			  lno_t,
			  scalar_t,
	          ExecSpace, ExecSpace, ExecSpace > KernelHandle;

	  KernelHandle kh;

	  kh.create_gs_handle(/*KokkosSparse::GS_TEAM*/);
	  Kokkos::Impl::Timer timer1;
	  KokkosKernels::Experimental::Example::pcgsolve(
	        kh
	      , crsmat
	      , kok_b_vector
	      , kok_x_vector
	      , cg_iteration_limit
	      , cg_iteration_tolerance
	      , & cg_result
	      , true
	  );
	  Kokkos::fence();

	  solve_time = timer1.seconds();


	  std::cout  << "DEFAULT SOLVE:"
	      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
	      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
	      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
	      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
	      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
	      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
	      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
	      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
	      << std::endl ;


#if KOKKOSSPARSE_IMPL_PRINTDEBUG
	  kok_x_vector = scalar_view_t("kok_x_vector", nv);

	  kh.create_gs_handle(KokkosSparse::GS_TEAM);
	  timer1.reset();
	  KokkosKernels::Experimental::Example::pcgsolve(
	        kh
	      , crsmat
	      , kok_b_vector
	      , kok_x_vector
	      , cg_iteration_limit
	      , cg_iteration_tolerance
	      , & cg_result
	      , true
	  );
	  Kokkos::fence();

	  solve_time = timer1.seconds();



	  std::cout  << "TEAM SOLVE:"
	      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
	      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
	      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
	      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
	      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
	      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
	      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
	      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
	      << std::endl ;
#endif
}


template <typename ExecSpace, typename crsMat_t>
void run_block_experiment(
    crsMat_t point_crsmat,
	crsMat_t block_crsmat, int block_size,
	typename crsMat_t::values_type::non_const_type kok_x_original){


	  //typedef typename crsMat_t::StaticCrsGraphType graph_t;
	  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
	  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
	  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

	  typedef typename lno_nnz_view_t::value_type lno_t;
	  typedef typename lno_view_t::value_type size_type;
	  typedef typename scalar_view_t::value_type scalar_t;
	  INDEX_TYPE nv = point_crsmat.numRows();
	  //scalar_view_t kok_x_original = create_x_vector<scalar_view_t>(nv, MAXVAL);

	  //KokkosKernels::Impl::print_1Dview(kok_x_original);
	  scalar_view_t kok_b_vector = create_y_vector(point_crsmat, kok_x_original);

	  //create X vector
	  scalar_view_t kok_x_vector("kok_x_vector", nv);


	  double solve_time = 0;
	  //const unsigned cg_iteration_limit = 10;
	  const double   cg_iteration_tolerance     = 1e-7 ;

	  KokkosKernels::Experimental::Example::CGSolveResult cg_result ;




	  typedef KokkosKernels::Experimental::KokkosKernelsHandle
	        < size_type,
			  lno_t,
			  scalar_t,
	          ExecSpace, ExecSpace, ExecSpace > KernelHandle;

	  KernelHandle kh;


	  kh.create_gs_handle();
	  Kokkos::Impl::Timer timer1;
	  KokkosKernels::Experimental::Example::block_pcgsolve(
	        kh
	      , point_crsmat
		  , block_crsmat, block_size
	      , kok_b_vector
	      , kok_x_vector
	      , cg_iteration_limit
	      , cg_iteration_tolerance
	      , & cg_result
	      , true
	  );
	  Kokkos::fence();

	  solve_time = timer1.seconds();


	  std::cout  << "DEFAULT SOLVE:"
	      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
	      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
	      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
	      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
	      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
	      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
	      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
	      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
	      << std::endl ;
}


template <typename ExecSpace, typename crsMat_t>
void run_experiment(
    crsMat_t crsmat,
	typename  crsMat_t::values_type::non_const_type kok_x_original, int block_size = 5 ){

	//run_point_experiment<ExecSpace, crsMat_t>(crsmat, kok_x_original);

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

  lno_view_t pf_rm;
  lno_nnz_view_t pf_e;
  scalar_view_t pf_v;
  size_t out_r, out_c;

  //typedef typename lno_nnz_view_t::value_type lno_t;
  //typedef typename lno_view_t::value_type size_type;
  //typedef typename scalar_view_t::value_type scalar_t;
  KokkosKernels::Impl::kk_create_blockcrs_formated_point_crsmatrix(
		  block_size , crsmat.numRows(), crsmat.numCols(),
		  crsmat.graph.row_map, crsmat.graph.entries, crsmat.values,
		  out_r, out_c,
		  pf_rm, pf_e, pf_v);


#if KOKKOSSPARSE_IMPL_PRINTDEBUG
  std::cout << "nr:" << crsmat.numRows() << " nc:" << crsmat.numCols() << std::endl;

  KokkosKernels::Impl::print_1Dview(crsmat.graph.row_map);
  KokkosKernels::Impl::print_1Dview(crsmat.graph.entries);
  KokkosKernels::Impl::print_1Dview(crsmat.values);
  std::cout << "onr:" << out_r << " oc:" << out_c<< std::endl;

  KokkosKernels::Impl::print_1Dview(pf_rm);
  KokkosKernels::Impl::print_1Dview(pf_e);
  KokkosKernels::Impl::print_1Dview(pf_v);
#endif

  graph_t static_graph2 (pf_e, pf_rm);
  crsMat_t crsmat2("CrsMatrix2", out_c, pf_v, static_graph2);
  run_point_experiment<ExecSpace, crsMat_t>(crsmat2, kok_x_original);


  lno_view_t bf_rm;
  lno_nnz_view_t bf_e;
  scalar_view_t bf_v;
  size_t but_r, but_c;

  KokkosKernels::Impl::kk_create_blockcrs_from_blockcrs_formatted_point_crs(
		  block_size , out_r, out_c,
		  pf_rm, pf_e, pf_v,
		  but_r, but_c,
		  bf_rm, bf_e, bf_v);

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
  KokkosKernels::Impl::print_1Dview(bf_rm);
  KokkosKernels::Impl::print_1Dview(bf_e);
  KokkosKernels::Impl::print_1Dview(bf_v);
#endif
  graph_t static_graph3 (bf_e, bf_rm);
  crsMat_t crsmat3("CrsMatrix3", but_c, bf_v, static_graph3);
  run_block_experiment<ExecSpace, crsMat_t>(crsmat2, crsmat3, block_size, kok_x_original);

}




enum { CMD_USE_THREADS = 0
     , CMD_USE_NUMA
     , CMD_USE_CORE_PER_NUMA
     , CMD_USE_CUDA
     , CMD_USE_OPENMP
     , CMD_USE_CUDA_DEV
     , CMD_BIN_MTX
     , CMD_ERROR
     , CMD_COUNT };

int main (int argc, char ** argv){


  int cmdline[ CMD_COUNT ] ;
  char *mtx_bin_file = NULL;
  int block_size = 5;
  struct Kokkos::InitArguments kargs;

  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;


  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--threads" ) ) {
      kargs.num_threads = cmdline[ CMD_USE_THREADS ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--openmp" ) ) {
       kargs.num_threads = cmdline[ CMD_USE_OPENMP ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "--cuda" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
    }
    else if ( 0 == strcasecmp( argv[i] , "--mtx" ) ) {
      mtx_bin_file = argv[++i];
    }

    else if ( 0 == strcasecmp( argv[i] , "--block_size" ) ) {
    	block_size = atoi( argv[++i] );
    }

    else if ( 0 == strcasecmp( argv[i] , "--iteration" ) ) {
    	cg_iteration_limit = atoi( argv[++i] );
    }
    else {
      cmdline[ CMD_ERROR ] = 1 ;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      std::cerr << "OPTIONS\n\t--threads [numThreads]\n\t--openmp [numThreads]\n\t--cuda\n\t--cuda-dev[DeviceIndex]\n\t--mtx[binary_mtx_file]" << std::endl;

      return 0;
    }
  }

  if (mtx_bin_file == NULL){
    std::cerr << "Provide a mtx binary file" << std::endl ;
    std::cerr << "OPTIONS\n\t--threads [numThreads]\n\t--openmp [numThreads]\n\t--cuda\n\t--cuda-dev[DeviceIndex]\n\t--mtx[binary_mtx_file]" << std::endl;
    return 0;
  }
  std::cout << "Running experiments with block size:" << block_size << std::endl;


  Kokkos::initialize(kargs);


#if defined( KOKKOS_ENABLE_THREADS )

    if ( cmdline[ CMD_USE_THREADS ] ) {
      INDEX_TYPE nv = 0, ne = 0;
      INDEX_TYPE *xadj, *adj;
      SCALAR_TYPE *ew;


      KokkosKernels::Impl::read_matrix<INDEX_TYPE,INDEX_TYPE, SCALAR_TYPE> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);
      Kokkos::Threads::print_configuration(std::cout);

      typedef Kokkos::Threads myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<SCALAR_TYPE, INDEX_TYPE, myExecSpace, void, SIZE_TYPE > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
      typedef typename graph_t::entries_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);

      KokkosKernels::Impl::copy_vector<SCALAR_TYPE * , values_view_t, myExecSpace>(ne, ew, values_view);
      KokkosKernels::Impl::copy_vector<INDEX_TYPE * , cols_view_t, myExecSpace>(ne, adj, columns_view);
      KokkosKernels::Impl::copy_vector<INDEX_TYPE * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);
      delete [] xadj;
      delete [] adj;
      delete [] ew;

      values_view_t kok_x_original = create_x_vector<values_view_t>(((nv /block_size) + 1) * block_size, MAXVAL);
      for (INDEX_TYPE i = nv; i < ((nv /block_size) + 1) * block_size; ++i){
    	  kok_x_original(i) = 0;
      }
      run_experiment<myExecSpace, crsMat_t>(crsmat, kok_x_original, block_size);
    }

#endif

#if defined( KOKKOS_ENABLE_OPENMP )

    if ( cmdline[ CMD_USE_OPENMP ] ) {
      INDEX_TYPE nv = 0, ne = 0;
      INDEX_TYPE *xadj, *adj;
      SCALAR_TYPE *ew;


      Kokkos::OpenMP::print_configuration(std::cout);

      KokkosKernels::Impl::read_matrix<INDEX_TYPE,INDEX_TYPE, SCALAR_TYPE> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);


      typedef Kokkos::OpenMP myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<SCALAR_TYPE, INDEX_TYPE, myExecSpace, void, SIZE_TYPE > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
      typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);

      KokkosKernels::Impl::copy_vector<SCALAR_TYPE * , values_view_t, myExecSpace>(ne, ew, values_view);
      KokkosKernels::Impl::copy_vector<INDEX_TYPE * , cols_view_t, myExecSpace>(ne, adj, columns_view);
      KokkosKernels::Impl::copy_vector<INDEX_TYPE * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);

      //crsMat_t crsmat("CrsMatrix", nv, nv, ne, ew, xadj, adj);
      delete [] xadj;
      delete [] adj;
      delete [] ew;


      values_view_t kok_x_original = create_x_vector<values_view_t>(((nv /block_size) + 1) * block_size, MAXVAL);
      for (INDEX_TYPE i = nv; i < ((nv /block_size) + 1) * block_size; ++i){
    	  kok_x_original(i) = 0;
      }
      run_experiment<myExecSpace, crsMat_t>(crsmat, kok_x_original, block_size);

    }

#endif

#if defined( KOKKOS_ENABLE_CUDA )
    if ( cmdline[ CMD_USE_CUDA ] ) {
      // Use the last device:
      INDEX_TYPE nv = 0, ne = 0;
      INDEX_TYPE *xadj, *adj;
      SCALAR_TYPE *ew;
      Kokkos::Cuda::print_configuration(std::cout);

      KokkosKernels::Impl::read_matrix<INDEX_TYPE,INDEX_TYPE, SCALAR_TYPE> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);


      typedef Kokkos::Cuda myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<SCALAR_TYPE, INDEX_TYPE, myExecSpace, void, SIZE_TYPE > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
      typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);


      {
        typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
        typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
        typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

        for (INDEX_TYPE i = 0; i <= nv; ++i){
          hr(i) = xadj[i];
        }

        for (INDEX_TYPE i = 0; i < ne; ++i){
          hc(i) = adj[i];
          hv(i) = ew[i];
        }
        Kokkos::deep_copy (rowmap_view , hr);
        Kokkos::deep_copy (columns_view , hc);
        Kokkos::deep_copy (values_view , hv);


      }
      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);

      delete [] xadj;
      delete [] adj;
      delete [] ew;

      values_view_t kok_x_original = create_x_vector<values_view_t>(((nv /block_size) + 1) * block_size, MAXVAL);
      run_experiment<myExecSpace, crsMat_t>(crsmat, kok_x_original, block_size);


    }

#endif
  Kokkos::finalize();


  return 0;
}
#else //defined(KOKKOSKERNELS_INST_DOUBLE) &&  defined(KOKKOSKERNELS_INST_OFFSET_INT) &&     defined(KOKKOSKERNELS_INST_ORDINAL_INT)

int main() {
  std::cerr << "PCG is configured with INT, INT, DOUBLE. Test is not compiled as these are not instantiated." << std::endl; 
}
#endif
