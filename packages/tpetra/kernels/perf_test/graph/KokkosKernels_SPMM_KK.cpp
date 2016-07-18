
#include <iostream>
#include "KokkosKernels_GraphHelpers.hpp"
#include "KokkosKernels_SPGEMM.hpp"
#include <Kokkos_Sparse_CrsMatrix.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_GraphColor.hpp"

typedef int idx;
typedef double wt;


#define TRANPOSEFIRST false
#define TRANPOSESECOND false

template <typename ExecSpace, typename crsMat_t>
crsMat_t run_experiment(
    crsMat_t crsmat, crsMat_t crsmat2, int algorithm, int repeat , int chunksize, int multi_color_scale, int shmemsize, int teamsize, int use_dynamic_scheduling);
template <typename v1>
struct compare{
  v1 f,s;
  compare (v1 f_ , v1 s_): f(f_), s(s_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_t &diff) const {

    if (f[i] - s[i] > 0.00001 || f[i] - s[i] < -0.00001) diff++;
  }

};


enum { CMD_USE_THREADS = 0
  , CMD_USE_NUMA
  , CMD_USE_CORE_PER_NUMA
  , CMD_USE_CUDA
  , CMD_USE_OPENMP
  , CMD_USE_CUDA_DEV
  , CMD_SPGEMM_ALGO
  , CMD_BIN_AMTX
  , CMD_BIN_RMTX
  , CMD_BIN_PMTX
  , CMD_MM_MODE
  , CMD_REPEAT
  , CMD_CHUNKSIZE
  , CMD_MULTICOLORSCALE
  , CMD_SHMEMSIZE
  , CMD_TEAMSIZE
  , CMD_DYNAMIC_SCHEDULE
  , CMD_ERROR
  , CMD_COUNT };

int main (int argc, char ** argv){



  int cmdline[ CMD_COUNT ] ;
  char *r_mtx_bin_file = NULL;
  char *a_mtx_bin_file = NULL;
  char *p_mtx_bin_file = NULL;

  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;
  cmdline[ CMD_REPEAT ] = 1;
  cmdline[ CMD_CHUNKSIZE ] = -1;
  cmdline[ CMD_MULTICOLORSCALE ] = 1;
  cmdline[ CMD_SHMEMSIZE ] = 16128;
  cmdline[ CMD_TEAMSIZE ] = -1;
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "threads" ) ) {
      cmdline[ CMD_USE_THREADS ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) {
      cmdline[ CMD_USE_OPENMP ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "repeat" ) ) {
      cmdline[ CMD_REPEAT ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "cores" ) ) {
      sscanf( argv[++i] , "%dx%d" ,
          cmdline + CMD_USE_NUMA ,
          cmdline + CMD_USE_CORE_PER_NUMA );
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
    }
    else if ( 0 == strcasecmp( argv[i] , "chunksize" ) ) {
      cmdline[ CMD_CHUNKSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "teamsize" ) ) {
      cmdline[ CMD_TEAMSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda-dev" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
      cmdline[ CMD_USE_CUDA_DEV ] = atoi( argv[++i] ) ;
    }

    else if ( 0 == strcasecmp( argv[i] , "mmmode" ) ) {
      cmdline[ CMD_MM_MODE ] = atoi( argv[++i] ) ;
    }

    else if ( 0 == strcasecmp( argv[i] , "mcscale" ) ) {
      cmdline[ CMD_MULTICOLORSCALE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "shmem" ) ) {
      cmdline[ CMD_SHMEMSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "amtx" ) ) {
      a_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "rmtx" ) ) {
      r_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "pmtx" ) ) {
      p_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "dynamic" ) ) {
      cmdline[ CMD_DYNAMIC_SCHEDULE ]  = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "algorithm" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "KK" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 0;
      }
      else if ( 0 == strcasecmp( argv[i] , "MKL" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 1;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSPARSE" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 2;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSP" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 3;
      }
      else if ( 0 == strcasecmp( argv[i] , "DKK" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 4;
      }
      else if ( 0 == strcasecmp( argv[i] , "CKK" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 5;
      }
      else if ( 0 == strcasecmp( argv[i] , "SKK" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 6;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKSPEED" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 8;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEM" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 7;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKCOLOR" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 9;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 10;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR2" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 11;
      }

      else {
        cmdline[ CMD_ERROR ] = 1 ;
        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
        return 0;
      }
    }
    else {
      cmdline[ CMD_ERROR ] = 1 ;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;

      return 0;
    }
  }

  if (a_mtx_bin_file == NULL){
    std::cerr << "Provide a mtx binary file" << std::endl ;
    std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
    return 0;
  }
  if (cmdline[ CMD_MM_MODE ] == 1 || cmdline[ CMD_MM_MODE ] == 2){
    if (r_mtx_bin_file == NULL){
      std::cerr << "Provide a r mtx binary file rmtx rmatrix_file" << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
      return 0;
    }
    if (p_mtx_bin_file == NULL){
      std::cerr << "Provide a p mtx binary file pmtx rmatrix_file" << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
      return 0;
    }
  }

  idx m = 0, nnzA = 0;
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

    if (cmdline[ CMD_SPGEMM_ALGO ] == 2 || cmdline[ CMD_SPGEMM_ALGO ] == 3){
      std::cerr << "CUSP and CUSPARSE cannot be run with PTHREADS" << std::endl ;
      return 0;
    }

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = m;
    idx ne = nnzA;
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

    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
    }
    else if (cmdline[ CMD_MM_MODE ] == 1){
      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);


        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);

        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
    }


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
    if (cmdline[ CMD_SPGEMM_ALGO ] == 2 || cmdline[ CMD_SPGEMM_ALGO ] == 3){
      std::cerr << "CUSP and CUSPARSE cannot be run with OPENMP" << std::endl ;
      return 0;
    }

    Kokkos::OpenMP::print_configuration(std::cout);

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = m;
    idx ne = nnzA;


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
    crsMat_t crsmat("CrsMatrix", m, values_view, static_graph);

    delete [] xadj;
    delete [] adj;
    delete [] ew;


    std::cout << "STARTUP MULTIPLYING A*A" << std::endl;
    run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat,  cmdline[ CMD_SPGEMM_ALGO ], cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
    std::cout << "STARTUP DONE  A*A\n\n\n\n\n" << std::endl;

    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat,  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
    }else if (cmdline[ CMD_MM_MODE ] == 1){
      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);


        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);

        /*
        char *file = "AP.mtx";
        KokkosKernels::Experimental::Graph::Utils::write_graph_bin(
            m, (idx) crsmat.graph.entries.dimension_0(),
            ( const idx *)crsmat.graph.row_map.ptr_on_device(),
            ( const idx *)crsmat.graph.entries.ptr_on_device(),
            ( const wt *)crsmat.values.ptr_on_device(),
            ( const char *)file);
         */
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);


        delete [] xadj;
        delete [] adj;
        delete [] ew;

        //cmdline[ CMD_SPGEMM_ALGO ] = 1;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }

    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);


        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 3){
      {
        std::cout << "MULTIPLYING R*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        crsmat = crsmat2;

        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);
        rowmap_view = row_map_view_t ("rowmap_view", m+1);
        columns_view = cols_view_t ("colsmap_view", nnzA);
        values_view = values_view_t ("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        static_graph = graph_t (columns_view, rowmap_view);
        crsmat2 = crsMat_t ("CrsMatrix2", ncols, values_view, static_graph);

        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);


        if (0)
        {

          int *c_xadj = (int *)crsmat.graph.row_map.ptr_on_device();
          int *c_adj = (int *)crsmat.graph.entries.ptr_on_device();
          double *c_ew = (double *) crsmat.values.ptr_on_device();
          int m = crsmat.graph.row_map.dimension_0() - 1;
          int nnzA = crsmat.values.dimension_0();
          KokkosKernels::Experimental::Graph::Utils::write_graph_bin<int, double>
          (m, nnzA, c_xadj, c_adj, c_ew, "result.bin");
        }

      }
    }


    myExecSpace::finalize();
  }

#endif

#if defined( KOKKOS_HAVE_CUDA )
  if ( cmdline[ CMD_USE_CUDA ] ) {
    // Use the last device:

    if (cmdline[ CMD_SPGEMM_ALGO ] == 1){
      std::cerr << "MKL cannot be run with CUDA" << std::endl ;
      return 0;
    }
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cmdline[ CMD_USE_CUDA_DEV ] ) );
    Kokkos::Cuda::print_configuration(std::cout);

    {
      //just warm up gpu.
    Kokkos::View<int *, Kokkos::Cuda> tmp_view("rowmap_view", 20000000);
    Kokkos::deep_copy(tmp_view, 2);
    }



    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = m;
    idx ne = nnzA;

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
      typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
      typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
      typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

      for (idx i = 0; i <= nv; ++i){
        hr(i) = xadj[i];
      }

      for (idx i = 0; i < nnzA; ++i){
        hc(i) = adj[i];
        hv(i) = ew[i];
      }
      Kokkos::deep_copy (rowmap_view , hr);
      Kokkos::deep_copy (columns_view , hc);
      Kokkos::deep_copy (values_view , hv);


    }
    graph_t static_graph (columns_view, rowmap_view);
    crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);


    //n = m;


    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<Kokkos::Cuda, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
    }
    else if (cmdline[ CMD_MM_MODE ] == 1){

      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);

      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }


        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);



        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }

        //KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);

      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);
        std::cout << 1 << std::endl;

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);


        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);
          std::cout << 2 << std::endl;
          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }
          std::cout << 3 << std::endl;
          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          std::cout << 4 << std::endl;
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
          std::cout << 5 << std::endl;
        }

        //KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ]);
      }
    }




    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }

#endif


  return 0;

}


template <typename ExecSpace, typename crsMat_t>
crsMat_t run_experiment(
    crsMat_t crsMat, crsMat_t crsMat2,
    int algorithm, int repeat, int chunk_size ,int multi_color_scale, int shmemsize, int team_size, int use_dynamic_scheduling){


  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      ExecSpace, ExecSpace,ExecSpace > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);

  if (use_dynamic_scheduling){
    kh.set_dynamic_scheduling(true);
  }

  const size_t m = crsMat.numRows();
  const size_t n = crsMat2.numRows();
  const size_t k = crsMat2.numCols();

  std::cout << "m:" << m << " n:" << n << " k:" << k << std::endl;
  if (n != crsMat.numCols()){
    std::cout << "crsMat.numCols():" << crsMat.numCols() << " crsMat2.numRows():" << crsMat2.numRows() << std::endl;
    exit(1);
  }

  switch (algorithm){
  case 0:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK1);
    break;
  case 1:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_MKL);
    break;
  case 2:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_CUSPARSE);
    break;
  case 3:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_CUSP);
    break;
  case 4:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK2);
    break;
  case 5:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK3);
    break;
  case 6:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK4);
    break;
  case 7:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY);
      break;
  case 8:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED);
      break;
  case 9:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR);
    break;

  case 10:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR);
    break;

  case 11:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2);
      break;

  default:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY);
    break;
  }

  kh.get_spgemm_handle()->set_multi_color_scale(multi_color_scale);
  for (int i = 0; i < repeat; ++i){
    row_mapC = lno_view_t ("");
    entriesC = lno_nnz_view_t ("");
    valuesC = scalar_view_t ("");

    Kokkos::Impl::Timer timer1;
    KokkosKernels::Experimental::Graph::spgemm_symbolic (
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        TRANPOSEFIRST,
        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        TRANPOSESECOND,
        row_mapC,
        entriesC
    );

    ExecSpace::fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Impl::Timer timer2;
    KokkosKernels::Experimental::Graph::spgemm_numeric(
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        crsMat.values,
        TRANPOSEFIRST,

        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        crsMat2.values,
        TRANPOSESECOND,
        row_mapC,
        entriesC,
        valuesC
    );
    ExecSpace::fence();
    double numeric_time = timer2.seconds();

    Kokkos::Impl::Timer timer3;
    KokkosKernels::Experimental::Graph::spgemm_apply(
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        crsMat.values,
        TRANPOSEFIRST,

        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        crsMat2.values,
        TRANPOSESECOND,
        row_mapC,
        entriesC,
        valuesC
    );
    ExecSpace::fence();
    double apply_time = timer3.seconds();

    std::cout
    << "mm_time:" << numeric_time + symbolic_time + apply_time
    << " symbolic_time:" << symbolic_time
    << " numeric:" << numeric_time
    << " apply:" << apply_time << std::endl;
  }

  std::cout << "row_mapC:" << row_mapC.dimension_0() << std::endl;
  std::cout << "entriesC:" << entriesC.dimension_0() << std::endl;
  std::cout << "valuesC:" << valuesC.dimension_0() << std::endl;



  {


    typename lno_view_t::HostMirror hr = Kokkos::create_mirror_view (row_mapC);
    Kokkos::deep_copy ( hr, row_mapC);

    ExecSpace::fence();


    idx *cent = entriesC.ptr_on_device();
    wt * cval = valuesC.ptr_on_device();

    scalar_view_t my_row_0_vals ("row0", hr(1));
    lno_nnz_view_t my_row_0_entries ("row0", hr(1));

    ExecSpace::fence();
    KokkosKernels::Experimental::Util::copy_vector<wt * , scalar_view_t, ExecSpace>(hr(1), cval, my_row_0_vals);
    ExecSpace::fence();
    KokkosKernels::Experimental::Util::copy_vector<idx * , lno_nnz_view_t, ExecSpace>(hr(1), cent, my_row_0_entries);
    ExecSpace::fence();

    KokkosKernels::Experimental::KokkosKernelsSorting::sort_key_value_views<lno_nnz_view_t,scalar_view_t, ExecSpace>(my_row_0_entries, my_row_0_vals);

    std::cout << "RESULT FIRST ROW" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(my_row_0_entries, true);
    KokkosKernels::Experimental::Util::print_1Dview(my_row_0_vals, true);
    std::cout << "#################" << std::endl;
  }

  if (0)
  {
    typename lno_view_t::HostMirror hr = Kokkos::create_mirror_view (row_mapC);
    Kokkos::deep_copy ( hr, row_mapC);

    typename lno_nnz_view_t::HostMirror he = Kokkos::create_mirror_view (entriesC);
    Kokkos::deep_copy ( he, entriesC);

    typename scalar_view_t::HostMirror hv = Kokkos::create_mirror_view (valuesC);
    Kokkos::deep_copy ( hv, valuesC);

    std::vector <KokkosKernels::Experimental::Graph::Utils::Edge<idx, wt>> edge_list (he.dimension_0());


    for (idx i = 0; i < hr.dimension_0() - 1; ++i){
      idx begin = hr(i);
      idx end = hr(i + 1);
      idx edge_ind = 0;
      for (idx j = begin; j < end; ++j){
        edge_list[edge_ind].src = i;
        edge_list[edge_ind].dst = he(j);
        edge_list[edge_ind].ew = hv(j);
        edge_ind++;
      }
      std::sort (edge_list.begin(), edge_list.begin() + edge_ind);

      edge_ind = 0;
      for (int j = begin; j < end; ++j){
        he(j) = edge_list[edge_ind].dst;
        hv(j) = edge_list[edge_ind++].ew;
      }

      Kokkos::deep_copy ( entriesC, he);
      Kokkos::deep_copy ( valuesC, hv);
    }

    KokkosKernels::Experimental::Util::print_1Dview(entriesC);
    KokkosKernels::Experimental::Util::print_1Dview(valuesC);
  }
/*
  if (1)
  {
    std::cout << "D1 Coloring Result Matrix " << std::endl;
    kh.create_graph_coloring_handle();
    typename KernelHandle::GraphColoringHandleType *gch = kh.get_graph_coloring_handle();
    gch->set_coloring_type(KokkosKernels::Experimental::Graph::Distance1);
    gch->set_algorithm(KokkosKernels::Experimental::Graph::COLORING_SERIAL);

    lno_view_t tmp_xadj;
    lno_nnz_view_t tmp_adj;


    Kokkos::Impl::Timer timer;
    KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap
    < lno_view_t, lno_nnz_view_t,
    lno_view_t, lno_nnz_view_t,
    ExecSpace>
    (m, row_mapC, entriesC, tmp_xadj, tmp_adj );



    std::cout << " Symmetrized Graph: NV:" << m << " NE:" << entriesC.dimension_0() << " symmetrization time:" << timer.seconds() << std::endl;
    timer.reset();

    KokkosKernels::Experimental::Graph::graph_color_symbolic
        <KernelHandle,lno_view_t,lno_nnz_view_t> (&kh,m, m , tmp_xadj, tmp_adj);




    std::cout << "Num colors:" << gch->get_num_colors() <<  " coloring time:" << timer.seconds()  << std::endl;


    typename KernelHandle::GraphColoringHandleType::color_view_t color_view = kh.get_graph_coloring_handle()->get_vertex_colors();

    lno_view_t histogram ("histogram", gch->get_num_colors());

    ExecSpace::fence();

    timer.reset();


    std::cout << "Histogram" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);
    std::cout << "Colors" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(color_view);

    KokkosKernels::Experimental::Util::get_histogram
      <typename KernelHandle::GraphColoringHandleType::color_view_t, lno_view_t, ExecSpace>(m, color_view, histogram);


    std::cout << "Histogram" << " time:" << timer.seconds()  << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);

    kh.destroy_graph_coloring_handle();
  }

  if (1)
  {
    std::cout << "Coloring Result Matrix with new distance-2 color" << std::endl;
    kh.create_graph_coloring_handle();
    typename KernelHandle::GraphColoringHandleType *gch = kh.get_graph_coloring_handle();
    gch->set_coloring_type(KokkosKernels::Experimental::Graph::Distance2);
    gch->set_algorithm(KokkosKernels::Experimental::Graph::COLORING_SERIAL2);

    lno_view_t tmp_xadj;
    lno_nnz_view_t tmp_adj;


    Kokkos::Impl::Timer timer;
    KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap
    < lno_view_t, lno_nnz_view_t,
    lno_view_t, lno_nnz_view_t,
    ExecSpace>
    (m, row_mapC, entriesC, tmp_xadj, tmp_adj );



    std::cout << " Symmetrized Graph: NV:" << m << " NE:" << entriesC.dimension_0() << " symmetrization time:" << timer.seconds() << std::endl;
    timer.reset();

    KokkosKernels::Experimental::Graph::graph_color_symbolic
        <KernelHandle,lno_view_t,lno_nnz_view_t> (&kh,m, m , tmp_xadj, tmp_adj);




    std::cout << "Num colors:" << gch->get_num_colors() <<  " coloring time:" << timer.seconds()  << std::endl;


    typename KernelHandle::GraphColoringHandleType::color_view_t color_view = kh.get_graph_coloring_handle()->get_vertex_colors();

    lno_view_t histogram ("histogram", gch->get_num_colors());

    ExecSpace::fence();

    timer.reset();


    std::cout << "Histogram" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);
    std::cout << "Colors" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(color_view);

    KokkosKernels::Experimental::Util::get_histogram
      <typename KernelHandle::GraphColoringHandleType::color_view_t, lno_view_t, ExecSpace>(m, color_view, histogram);


    std::cout << "Histogram" << " time:" << timer.seconds()  << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);

    kh.destroy_graph_coloring_handle();
  }



  if (1)
  {
    std::cout << "Coloring Result Matrix" << std::endl;
    kh.create_graph_coloring_handle();
    typename KernelHandle::GraphColoringHandleType *gch = kh.get_graph_coloring_handle();
    gch->set_coloring_type(KokkosKernels::Experimental::Graph::Distance2);

    lno_view_t tmp_xadj;
    lno_nnz_view_t tmp_adj;


    Kokkos::Impl::Timer timer;
    KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap
    < lno_view_t, lno_nnz_view_t,
    lno_view_t, lno_nnz_view_t,
    ExecSpace>
    (m, row_mapC, entriesC, tmp_xadj, tmp_adj );



    std::cout << " Symmetrized Graph: NV:" << m << " NE:" << entriesC.dimension_0() << " symmetrization time:" << timer.seconds() << std::endl;
    timer.reset();

    KokkosKernels::Experimental::Graph::graph_color_symbolic
        <KernelHandle,lno_view_t,lno_nnz_view_t> (&kh,m, m , tmp_xadj, tmp_adj);




    std::cout << "Num colors:" << gch->get_num_colors() <<  " coloring time:" << timer.seconds()  << std::endl;


    typename KernelHandle::GraphColoringHandleType::color_view_t color_view = kh.get_graph_coloring_handle()->get_vertex_colors();

    lno_view_t histogram ("histogram", gch->get_num_colors());

    ExecSpace::fence();

    timer.reset();


    std::cout << "Histogram" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);
    std::cout << "Colors" << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(color_view);

    KokkosKernels::Experimental::Util::get_histogram
      <typename KernelHandle::GraphColoringHandleType::color_view_t, lno_view_t, ExecSpace>(m, color_view, histogram);


    std::cout << "Histogram" << " time:" << timer.seconds()  << std::endl;
    KokkosKernels::Experimental::Util::print_1Dview(histogram);

    kh.destroy_graph_coloring_handle();
  }
*/


  typename crsMat_t::StaticCrsGraphType static_graph (entriesC, row_mapC);
  crsMat_t Ccrsmat("CrsMatrixC", k, valuesC, static_graph);
  return Ccrsmat;

}


