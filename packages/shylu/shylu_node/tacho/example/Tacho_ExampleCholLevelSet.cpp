#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Internal.hpp"
#include "Tacho_CommandLineParser.hpp"

//#define TACHO_ENABLE_PROFILE
#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
#include "cuda_profiler_api.h" 
#endif

//#define TACHO_ENABLE_MPI_TEST
#if defined(TACHO_ENABLE_MPI_TEST)
#include "mpi.h"
#endif

template<typename T> using TaskSchedulerType = Kokkos::TaskSchedulerMultiple<T>;
static const char * scheduler_name = "TaskSchedulerMultiple";

int main (int argc, char *argv[]) {
  Tacho::CommandLineParser opts("This example program measure the performance of Tacho level set tools");

  bool serial = false;
  int nthreads = 1;
  bool verbose = true;
  bool sanitize = false;
  std::string file = "test.mtx";
  int nrhs = 1;
  int device_level_cut = 0;
  int device_factor_thres = 256;
  int device_solve_thres = 256;
  int nstreams = 8;

  opts.set_option<bool>("serial", "Flag to use serial algorithm", &serial);
  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("nstreams", "# of streams used in level set factorization and solve", &nstreams);
  opts.set_option<int>("device-level-cut", "tree cut level to force device function", &device_level_cut);
  opts.set_option<int>("device-factor-thres", "device function is used above this threshold in factorization", &device_factor_thres);
  opts.set_option<int>("device-solve-thres", "device function is used above this threshold in solve", &device_solve_thres);

#if defined(TACHO_ENABLE_MPI_TEST)
  int nrank, irank;
  MPI_Init(&argc, &argv);
  const auto comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &irank);
  MPI_Comm_size(comm, &nrank);

  std::vector<double> t_all(nrank, double(0));
  auto is_root = [irank]()->bool { return irank == 0; };
#else
  auto is_root = []()->bool { return true; };
#endif

#if !defined (KOKKOS_ENABLE_CUDA)
  // override serial flag
  if (serial) {
    std::cout << "CUDA is enabled and serial code cannot be instanciated\n";
    serial = false;
  }
#endif

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);

#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
  cudaProfilerStop();
#endif

  typedef typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::device_type device_type;
  typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;

  typedef TaskSchedulerType<typename device_type::execution_space> scheduler_type;

  Tacho::printExecSpaceConfiguration<typename device_type::execution_space> ("DeviceSpace", false);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space> ("HostSpace", false);
  printf("Scheduler Type = %s\n", scheduler_name);
  
  int r_val = 0;
  
  {
    typedef double value_type;
    typedef Tacho::CrsMatrixBase<value_type,device_type> CrsMatrixBaseType;
    typedef Tacho::CrsMatrixBase<value_type,host_device_type> CrsMatrixBaseTypeHost;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> DenseMatrixBaseType;
    
    Kokkos::Impl::Timer timer;
    double t = 0.0;
    
    if (is_root()) std::cout << "CholLevelSet:: import input file = " << file << std::endl;
    CrsMatrixBaseTypeHost A;
    timer.reset();
    {
      {
        std::ifstream in;
        in.open(file);
        if (!in.good()) {
          std::cout << "Failed in open the file: " << file << std::endl;
          return -1;
        }
      }
      Tacho::MatrixMarket<value_type>::read(file, A, sanitize, verbose);
    }
    Tacho::Graph G(A);
    t = timer.seconds();
    if (is_root()) std::cout << "CholLevelSet:: import input file::time = " << t << std::endl;

    if (is_root()) std::cout << "CholLevelSet:: analyze matrix" << std::endl;
    timer.reset();
#if   defined(TACHO_HAVE_METIS)
    Tacho::GraphTools_Metis T(G);
#elif defined(TACHO_HAVE_SCOTCH)
    Tacho::GraphTools_Scotch T(G);
#else
    Tacho::GraphTools_CAMD T(G);
#endif
    T.reorder(verbose);
    
    Tacho::SymbolicTools S(A, T);
    S.symbolicFactorize(verbose);
    t = timer.seconds();
    if (is_root()) std::cout << "CholLevelSet:: analyze matrix::time = " << t << std::endl;

    typedef typename device_type::memory_space device_memory_space; 
    
    auto a_row_ptr              = Kokkos::create_mirror_view(device_memory_space(), A.RowPtr());
    auto a_cols                 = Kokkos::create_mirror_view(device_memory_space(), A.Cols());
    auto a_values               = Kokkos::create_mirror_view(device_memory_space(), A.Values());
    
    auto t_perm                 = Kokkos::create_mirror_view(device_memory_space(), T.PermVector());
    auto t_peri                 = Kokkos::create_mirror_view(device_memory_space(), T.InvPermVector());
    auto s_supernodes           = Kokkos::create_mirror_view(device_memory_space(), S.Supernodes());
    auto s_gid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelPtr());
    auto s_gid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelColIdx());
    auto s_sid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelPtr());
    auto s_sid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelColIdx());
    auto s_blk_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.blkSuperPanelColIdx());
    auto s_snodes_tree_parent   = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeParent());
    auto s_snodes_tree_ptr      = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreePtr());
    auto s_snodes_tree_children = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeChildren());
    
    Kokkos::deep_copy(a_row_ptr              , A.RowPtr());
    Kokkos::deep_copy(a_cols                 , A.Cols());
    Kokkos::deep_copy(a_values               , A.Values());
    
    Kokkos::deep_copy(t_perm                 , T.PermVector());
    Kokkos::deep_copy(t_peri                 , T.InvPermVector());
    Kokkos::deep_copy(s_supernodes           , S.Supernodes());
    Kokkos::deep_copy(s_gid_spanel_ptr       , S.gidSuperPanelPtr());
    Kokkos::deep_copy(s_gid_spanel_colidx    , S.gidSuperPanelColIdx());
    Kokkos::deep_copy(s_sid_spanel_ptr       , S.sidSuperPanelPtr());
    Kokkos::deep_copy(s_sid_spanel_colidx    , S.sidSuperPanelColIdx());
    Kokkos::deep_copy(s_blk_spanel_colidx    , S.blkSuperPanelColIdx());
    Kokkos::deep_copy(s_snodes_tree_parent   , S.SupernodesTreeParent());
    Kokkos::deep_copy(s_snodes_tree_ptr      , S.SupernodesTreePtr());
    Kokkos::deep_copy(s_snodes_tree_children , S.SupernodesTreeChildren());

    Tacho::NumericTools<value_type,scheduler_type> 
      N(A.NumRows(), a_row_ptr, a_cols,
        t_perm, t_peri,
        S.NumSupernodes(), s_supernodes,
        s_gid_spanel_ptr, s_gid_spanel_colidx,
        s_sid_spanel_ptr, s_sid_spanel_colidx, s_blk_spanel_colidx,
        s_snodes_tree_parent, s_snodes_tree_ptr, s_snodes_tree_children, 
        S.SupernodesTreeLevel(),
        S.SupernodesTreeRoots());
    N.printMemoryStat(verbose);

    Tacho::LevelSetTools<value_type,scheduler_type> L(N);
    L.initialize(device_level_cut, device_factor_thres, device_solve_thres, verbose);
    L.createStream(nstreams);


    if (is_root()) std::cout << "CholLevelSet:: factorize matrix" << std::endl;
#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
    cudaProfilerStart();
#endif
    timer.reset();    
    L.factorizeCholesky(a_values, verbose);
    t = timer.seconds();    
#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
    cudaProfilerStop();
#endif
    if (is_root()) std::cout << "CholLevelSet:: factorize matrix::time = " << t << std::endl;

#if defined(TACHO_ENABLE_MPI_TEST)
    {
      MPI_Gather(&t, 1, MPI_DOUBLE, t_all.data(), 1, MPI_DOUBLE, 0, comm); 
      
      if (is_root()) {
        double t_min(1000000), t_max(0), t_avg(0);
        for (int i=0;i<nrank;++i) {
          t_min = std::min(t_min, t_all[i]);
          t_max = std::max(t_max, t_all[i]);
          t_avg += t_all[i];
        }
        t_avg /= double(nrank);
        std::cout << "CholLevelSet:: factorize matrix::time (min, max, avg) = " << t_min << ", " << t_max << ", " << t_avg << "\n";  
      }
    }
#endif

    DenseMatrixBaseType 
      B("B", A.NumRows(), nrhs), 
      X("X", A.NumRows(), nrhs), 
      W("W", A.NumRows(), nrhs);

    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);      
      Kokkos::fill_random(B, random, value_type(1));
    }

    ///
    /// solve via level set
    ///
    if (is_root()) std::cout << "CholLevelSet:: solve matrix via LevelSetTools" << std::endl;

    timer.reset();        
#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
    cudaProfilerStart();
#endif
    constexpr int niter = 10;
    for (int i=0;i<niter;++i) 
      L.solveCholesky(X, B, W, verbose); 
#if defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_ENABLE_PROFILE)
    cudaProfilerStop();
#endif
    t = timer.seconds();    
    if (is_root()) std::cout << "CholLevelSet:: solve matrix::time = " << t << std::endl;
    L.release(verbose);

#if defined(TACHO_ENABLE_MPI_TEST)
    {
      MPI_Gather(&t, 1, MPI_DOUBLE, t_all.data(), 1, MPI_DOUBLE, 0, comm); 
      
      if (is_root()) {
        double t_min(1000000), t_max(0), t_avg(0);
        for (int i=0;i<nrank;++i) {
          t_min = std::min(t_min, t_all[i]);
          t_max = std::max(t_max, t_all[i]);
          t_avg += t_all[i];
        }
        t_avg /= double(nrank*niter);
        t_min /= double(niter);
        t_max /= double(niter);
        std::cout << "CholLevelSet:: solve matrix::time (min, max, avg) = " << t_min << ", " << t_max << ", " << t_avg << "\n";  
      }
    }
#endif

    CrsMatrixBaseType AA;
    AA.createMirror(A);
    AA.copy(A);

    const double res = Tacho::computeRelativeResidual(AA, X, B);
    if (is_root()) std::cout << "CholLevelSet:: residual = " << res << std::endl;

    N.release(verbose);
  }
  Kokkos::finalize();

#if defined(TACHO_ENABLE_MPI_TEST)
  MPI_Finalize();
#endif

  return r_val;
}
