#include "KokkosKernels_PCG.hpp"
#include "KokkosKernels_GraphHelpers.hpp"
#include "Kokkos_Sparse_MV.hpp"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include "KokkosKernels_Utils.hpp"
#include <iostream>

#define MAXVAL 1

typedef int idx;
typedef double wt;



template<typename scalar_view_t>
scalar_view_t create_x_vector(idx nv, wt max_value = 1.0){
  scalar_view_t kok_x ("X", nv);

  typename scalar_view_t::HostMirror h_x =  Kokkos::create_mirror_view (kok_x);


  for (idx i = 0; i < nv; ++i){
    wt r = static_cast <wt> (rand()) / static_cast <wt> (RAND_MAX / max_value);
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
void run_experiment(
    crsMat_t crsmat){

  typedef typename crsMat_t::values_type scalar_view_t;

  idx nv = crsmat.numRows();
  scalar_view_t kok_x_original = create_x_vector<scalar_view_t>(nv, MAXVAL);
  scalar_view_t kok_b_vector = create_y_vector(crsmat, kok_x_original);

  //create X vector
  scalar_view_t kok_x_vector("kok_x_vector", nv);


  double solve_time = 0;
  const unsigned cg_iteration_limit = 100000;
  const double   cg_iteration_tolerance     = 1e-7 ;

  KokkosKernels::Experimental::Example::CGSolveResult cg_result ;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
        < typename crsMat_t::StaticCrsGraphType::row_map_type,
          typename crsMat_t::StaticCrsGraphType::entries_type,
          typename crsMat_t::values_type,
          ExecSpace, ExecSpace, ExecSpace > KernelHandle;

  KernelHandle kh;


  kh.create_gs_handle();
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



  /*
  kh.destroy_gs_handle();
  kh.create_gs_handle(KokkosKernels::Experimental::Graph::GS_PERMUTED);

  kok_x_vector = scalar_view_t("kok_x_vector", nv);
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
  std::cout  << "\nPERMUTED SGS SOLVE:"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;


  kh.destroy_gs_handle();
  kh.create_gs_handle(KokkosKernels::Experimental::Graph::GS_TEAM);

  kok_x_vector = scalar_view_t("kok_x_vector", nv);
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
  std::cout  << "\nTEAM SGS SOLVE:"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;




  kok_x_vector = scalar_view_t("kok_x_vector", nv);
  timer1.reset();
  KokkosKernels::Experimental::Example::pcgsolve(
        kh
      , crsmat
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , false
  );
  Kokkos::fence();

  solve_time = timer1.seconds();
  std::cout  << "\nCG SOLVE (With no Preconditioner):"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;
  */
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
  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;


  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "threads" ) ) {
      cmdline[ CMD_USE_THREADS ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) {
      cmdline[ CMD_USE_OPENMP ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "cores" ) ) {
      sscanf( argv[++i] , "%dx%d" ,
              cmdline + CMD_USE_NUMA ,
              cmdline + CMD_USE_CORE_PER_NUMA );
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda-dev" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
      cmdline[ CMD_USE_CUDA_DEV ] = atoi( argv[++i] ) ;
    }

    else if ( 0 == strcasecmp( argv[i] , "mtx" ) ) {
      mtx_bin_file = argv[++i];
    }
    else {
      cmdline[ CMD_ERROR ] = 1 ;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[mtx][binary_mtx_file]" << std::endl;

      return 0;
    }
  }

  if (mtx_bin_file == NULL){
    std::cerr << "Provide a mtx binary file" << std::endl ;
    std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[mtx][binary_mtx_file]" << std::endl;

    return 0;
  }

  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;


#if defined( KOKKOS_HAVE_PTHREAD )

    if ( cmdline[ CMD_USE_THREADS ] ) {

      if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
        Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] ,
                                     cmdline[ CMD_USE_NUMA ] ,
                                     cmdline[ CMD_USE_CORE_PER_NUMA ] );
      }
      else {
        Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] );
      }

      KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);
      Kokkos::Threads::print_configuration(std::cout);

      typedef Kokkos::Threads myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, idx > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
      typedef typename graph_t::entries_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);

      KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(ne, ew, values_view);
      KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(ne, adj, columns_view);
      KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);
      delete [] xadj;
      delete [] adj;
      delete [] ew;

      run_experiment<myExecSpace, crsMat_t>(crsmat);

      myExecSpace::finalize();
    }

#endif

#if defined( KOKKOS_HAVE_OPENMP )

    if ( cmdline[ CMD_USE_OPENMP ] ) {

      if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
        Kokkos::OpenMP::initialize( cmdline[ CMD_USE_OPENMP ] ,
                                     cmdline[ CMD_USE_NUMA ] ,
                                     cmdline[ CMD_USE_CORE_PER_NUMA ] );
      }
      else {
        Kokkos::OpenMP::initialize( cmdline[ CMD_USE_OPENMP ] );
      }
      Kokkos::OpenMP::print_configuration(std::cout);

      KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);


      typedef Kokkos::OpenMP myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, idx > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
      typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);

      KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(ne, ew, values_view);
      KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(ne, adj, columns_view);
      KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);

      //crsMat_t crsmat("CrsMatrix", nv, nv, ne, ew, xadj, adj);
      delete [] xadj;
      delete [] adj;
      delete [] ew;

      run_experiment<myExecSpace, crsMat_t>(crsmat);

      myExecSpace::finalize();
    }

#endif

#if defined( KOKKOS_HAVE_CUDA )
    if ( cmdline[ CMD_USE_CUDA ] ) {
      // Use the last device:

      Kokkos::HostSpace::execution_space::initialize();
      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cmdline[ CMD_USE_CUDA_DEV ] ) );
      Kokkos::Cuda::print_configuration(std::cout);

      KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&nv, &ne, &xadj, &adj, &ew, mtx_bin_file);


      typedef Kokkos::Cuda myExecSpace;
      typedef typename KokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, idx > crsMat_t;

      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
      typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
      typedef typename crsMat_t::values_type::non_const_type values_view_t;

      row_map_view_t rowmap_view("rowmap_view", nv+1);
      cols_view_t columns_view("colsmap_view", ne);
      values_view_t values_view("values_view", ne);


      {
        typename row_map_view_t::rowmap_view hr = Kokkos::create_mirror_view (rowmap_view);
        typename cols_view_t::rowmap_view hc = Kokkos::create_mirror_view (columns_view);
        typename values_view_t::rowmap_view hv = Kokkos::create_mirror_view (values_view);

        for (idx i = 0; i <= nv; ++i){
          hr(i) = xadj[i];
        }

        for (idx i = 0; i < ne; ++i){
          hc(i) = adj[i];
          hv(i) = ew[i];
        }
        Kokkos::deep_copy (rowmap_view , hr);
        Kokkos::deep_copy (columns_view , hc);
        Kokkos::deep_copy (values_view , hv);


      }
      graph_t static_graph (columns_view, rowmap_view);
      crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);

//      typedef typename KokkosSparse::CrsMatrix<wt, idx, Kokkos::Cuda> crsMat_t;
//      crsMat_t crsmat("CrsMatrix", nv, nv, ne, ew, xadj, adj);
      delete [] xadj;
      delete [] adj;
      delete [] ew;


      run_experiment<myExecSpace, crsMat_t>(crsmat);

      myExecSpace::finalize();
      Kokkos::HostSpace::execution_space::finalize();
    }

#endif



  return 0;
}
