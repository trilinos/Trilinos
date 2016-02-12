
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
void run_experiment(
    crsMat_t crsmat);
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

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, mtx_bin_file);

    n = m;
    Kokkos::Threads::print_configuration(std::cout);
    typedef typename KokkosSparse::CrsMatrix<wt, idx, Kokkos::Threads> crsMat_t;
    crsMat_t crsmat("CrsMatrix", m, n, nnzA, ew, xadj, adj);
    delete [] xadj;
    delete [] adj;
    delete [] ew;


    run_experiment<Kokkos::Threads, crsMat_t>(crsmat);

    Kokkos::Threads::finalize();
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

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, mtx_bin_file);

    n = m;
    typedef typename KokkosSparse::CrsMatrix<wt, idx, Kokkos::OpenMP> crsMat_t;
    crsMat_t crsmat("CrsMatrix", m, n, nnzA, ew, xadj, adj);
    delete [] xadj;
    delete [] adj;
    delete [] ew;

    run_experiment<Kokkos::OpenMP, crsMat_t>(crsmat);

    Kokkos::OpenMP::finalize();
  }

#endif

#if defined( KOKKOS_HAVE_CUDA )
  if ( cmdline[ CMD_USE_CUDA ] ) {
    // Use the last device:

    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cmdline[ CMD_USE_CUDA_DEV ] ) );
    Kokkos::Cuda::print_configuration(std::cout);

    KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, mtx_bin_file);

    n = m;
    typedef typename KokkosSparse::CrsMatrix<wt, idx, Kokkos::Cuda> crsMat_t;
    crsMat_t crsmat("CrsMatrix", m, n, nnzA, ew, xadj, adj);
    delete [] xadj;
    delete [] adj;
    delete [] ew;


    run_experiment<Kokkos::Cuda, crsMat_t>(crsmat);

    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }

#endif


  return 0;

}


template <typename ExecSpace, typename crsMat_t>
void run_experiment(
    crsMat_t crsMat){

  typedef typename crsMat_t::values_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type lno_nnz_view_t;

  lno_view_t row_mapC, row_mapC2;
  lno_nnz_view_t entriesC, entriesC2;
  scalar_view_t valuesC, valuesC2;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      ExecSpace, ExecSpace,ExecSpace > KernelHandle;

  KernelHandle kh;


  const size_t m = crsMat.numRows();

  size_t n = m, k = m;


  kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK1);
  Kokkos::Impl::Timer timer1;
  KokkosKernels::Experimental::Graph::spgemm_symbolic (
      &kh,
      m,
      n,
      k,
      crsMat.graph.row_map,
      crsMat.graph.entries,
      TRANPOSEFIRST,
      crsMat.graph.row_map,
      crsMat.graph.entries,
      TRANPOSESECOND,
      row_mapC,
      entriesC
  );

  Kokkos::fence();
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

      crsMat.graph.row_map,
      crsMat.graph.entries,
      crsMat.values,
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

      crsMat.graph.row_map,
      crsMat.graph.entries,
      crsMat.values,
      TRANPOSESECOND,
      row_mapC,
      entriesC,
      valuesC
  );
  Kokkos::fence();
  double apply_time = timer3.seconds();

  std::cout << "mm_time:" << numeric_time + symbolic_time + apply_time
      << " symbolic_time:" << symbolic_time
      << " numeric:" << numeric_time
      << " apply:" << apply_time << std::endl;

  std::cout << "row_mapC:" << row_mapC.dimension_0() << std::endl;
  std::cout << "entriesC:" << entriesC.dimension_0() << std::endl;
  std::cout << "valuesC:" << valuesC.dimension_0() << std::endl;
}


