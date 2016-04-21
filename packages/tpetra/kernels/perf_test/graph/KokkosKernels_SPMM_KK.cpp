
#include <iostream>
#include "KokkosKernels_GraphHelpers.hpp"
#include "KokkosKernels_SPGEMM.hpp"
#include <Kokkos_Sparse_CrsMatrix.hpp>
#include "KokkosKernels_Handle.hpp"


typedef int idx;
typedef double wt;


#define TRANPOSEFIRST false
#define TRANPOSESECOND false

template <typename ExecSpace, typename crsMat_t>
crsMat_t run_experiment(
    crsMat_t crsmat, crsMat_t crsmat2, int algorithm);
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
     , CMD_ERROR
     , CMD_COUNT };

int main (int argc, char ** argv){



  int cmdline[ CMD_COUNT ] ;
  char *r_mtx_bin_file = NULL;
  char *a_mtx_bin_file = NULL;
  char *p_mtx_bin_file = NULL;

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

    else if ( 0 == strcasecmp( argv[i] , "mmmode" ) ) {
      cmdline[ CMD_MM_MODE ] = atoi( argv[++i] ) ;
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
      else {
        cmdline[ CMD_ERROR ] = 1 ;
        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
        return 0;
      }
    }
    else {
      cmdline[ CMD_ERROR ] = 1 ;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;

      return 0;
    }
  }

  if (a_mtx_bin_file == NULL){
    std::cerr << "Provide a mtx binary file" << std::endl ;
    std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
    return 0;
  }
  if (cmdline[ CMD_MM_MODE ] == 1 || cmdline[ CMD_MM_MODE ] == 2){
    if (r_mtx_bin_file == NULL){
      std::cerr << "Provide a r mtx binary file rmtx rmatrix_file" << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
      return 0;
    }
    if (p_mtx_bin_file == NULL){
      std::cerr << "Provide a p mtx binary file pmtx rmatrix_file" << std::endl ;
      std::cerr << "OPTIONS\n\tthreads [numThreads]\n\topenmp [numThreads]\n\tcuda\n\tcuda-dev[DeviceIndex]\n\t[[a|r|p]mtx][binary_mtx_file]\n\talgorithm[KK|MKL|CUSPARSE|CUSP]\n\t[mmode][0|1|2]: 0:AxA 1:R(AP) 2:(RA)P" << std::endl;
      return 0;
    }
  }

  idx m = 0, nnzA = 0, n = 0;
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
    idx nv = n = m;
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
      run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);

    }
    else if (cmdline[ CMD_MM_MODE ] == 1){
      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
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
    idx nv = n = m;
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

    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat,  cmdline[ CMD_SPGEMM_ALGO ]);
    }else if (cmdline[ CMD_MM_MODE ] == 1){
      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }

    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 3){
      {
        std::cout << "MULTIPLYING R*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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

        std::cout << "MULTIPLYING (RA)*P" << std::endl;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
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

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = n = m;
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


    n = m;


    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<Kokkos::Cuda, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
    }
    else if (cmdline[ CMD_MM_MODE ] == 1){

      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0, n = 0;
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
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ]);
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
    int algorithm){

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
  default:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK1);
    break;
  }
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

  Kokkos::fence();
  double symbolic_time = timer1.seconds();
  std::cout << "symbolic_time:" << symbolic_time << std::endl;
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
  Kokkos::fence();
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
  Kokkos::fence();
  double apply_time = timer3.seconds();

  std::cout
      << "mm_time:" << numeric_time + symbolic_time + apply_time
      << " symbolic_time:" << symbolic_time
      << " numeric:" << numeric_time
      << " apply:" << apply_time << std::endl;

  std::cout << "row_mapC:" << row_mapC.dimension_0() << std::endl;
  std::cout << "entriesC:" << entriesC.dimension_0() << std::endl;
  std::cout << "valuesC:" << valuesC.dimension_0() << std::endl;



  {
    for (int i = 0; i < row_mapC.dimension_0() - 1; ++i){
      auto subentries = Kokkos::subview(entriesC, Kokkos::make_pair(row_mapC(i), row_mapC(i + 1)));
      Kokkos::sort(subentries);
      //KokkosKernels::Experimental::Util::print_1Dview(subentries);
    }
  }


  typename crsMat_t::StaticCrsGraphType static_graph (entriesC, row_mapC);
  crsMat_t Ccrsmat("CrsMatrixC", k, valuesC, static_graph);
  return Ccrsmat;

}


