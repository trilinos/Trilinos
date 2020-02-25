#include "ShyLU_NodeTacho_config.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_Graph.hpp"
#include "Tacho_SymbolicTools.hpp"

#if defined(TACHO_HAVE_SCOTCH)
#include "Tacho_GraphTools_Scotch.hpp"
#endif

#if defined(TACHO_HAVE_METIS)
#include "Tacho_GraphTools_Metis.hpp"
#endif

#include "Tacho_GraphTools_CAMD.hpp"

#include "Tacho_NumericTools.hpp"
#include "Tacho_TriSolveTools.hpp"

#include "Tacho_CommandLineParser.hpp"

template<typename T> using TaskSchedulerType = Kokkos::TaskSchedulerMultiple<T>;
static const char * scheduler_name = "TaskSchedulerMultiple";

int main (int argc, char *argv[]) {
  Tacho::CommandLineParser opts("This example program measure the performance of Tacho on Kokkos::OpenMP");

  bool serial = false;
  int nthreads = 1;
  bool verbose = true;
  std::string file = "test.mtx";
  int max_num_superblocks = 4;
  int nrhs = 1;
  int serial_thres_size = -1; // 32 is better
  int mb = 0;
  int nb = 0;
  int device_level_cut = 0;
  int device_function_thres = 256;
  int nstreams_solve = 8;
  int nstreams_prepare = 2;

  opts.set_option<bool>("serial", "Flag to use serial algorithm", &serial);
  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("max-num-superblocks", "Max number of superblocks", &max_num_superblocks);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("serial-thres", "Serialization threshold size", &serial_thres_size);
  opts.set_option<int>("mb", "Blocksize", &mb);
  opts.set_option<int>("nb", "Panelsize", &nb);
  opts.set_option<int>("nstreams-solve", "# of streams used in trisolve", &nstreams_solve);
  opts.set_option<int>("nstreams-prepare", "# of streams used in trisolve", &nstreams_prepare);
  opts.set_option<int>("device-level-cut", "tree cut level to force device function", &device_level_cut);
  opts.set_option<int>("device-function-thres", "device function is used above this threshold", &device_function_thres);

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

  typedef typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::device_type device_type;
  typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;

  typedef TaskSchedulerType<typename device_type::execution_space> scheduler_type;

  Tacho::printExecSpaceConfiguration<typename device_type::execution_space> ("DeviceSpace", false);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space> ("HostSpace", false);
  printf("Scheduler Type = %s\n", scheduler_name);
  
  int r_val = 0;
  
  {
    typedef double value_type;
    typedef Tacho::CrsMatrixBase<value_type,host_device_type> CrsMatrixBaseTypeHost;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> DenseMatrixBaseType;
    
    Kokkos::Impl::Timer timer;
    double t = 0.0;

    std::cout << "CholTriSolve:: import input file = " << file << std::endl;
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
      Tacho::MatrixMarket<value_type>::read(file, A);
    }
    Tacho::Graph G(A);
    t = timer.seconds();
    std::cout << "CholTriSolve:: import input file::time = " << t << std::endl;

    std::cout << "CholTriSolve:: analyze matrix" << std::endl;
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
    std::cout << "CholTriSolve:: analyze matrix::time = " << t << std::endl;

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

    N.setSerialThresholdSize(serial_thres_size);
    N.setMaxNumberOfSuperblocks(max_num_superblocks);

    std::cout << "CholTriSolve:: factorize matrix" << std::endl;
    timer.reset();    
    if (serial) {
#if !defined (KOKKOS_ENABLE_CUDA)
      N.factorizeCholesky_Serial(a_values, verbose);
#endif
    } else {
      if (nb <= 0) {
        if (mb > 0)
          N.factorizeCholesky_ParallelByBlocks(a_values, mb, verbose);
        else
          N.factorizeCholesky_Parallel(a_values, verbose);
      } else {
        if (mb > 0)
          N.factorizeCholesky_ParallelByBlocksPanel(a_values, mb, nb, verbose);
        else
          N.factorizeCholesky_ParallelPanel(a_values, nb, verbose);
      }
    }
    t = timer.seconds();    
    std::cout << "CholTriSolve:: factorize matrix::time = " << t << std::endl;
    
    DenseMatrixBaseType 
      B("B", A.NumRows(), nrhs), 
      X("X", A.NumRows(), nrhs), 
      W("W", A.NumRows(), nrhs);

    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);      
      Kokkos::fill_random(B, random, value_type(1));
    }

#if 1
    ///
    /// solve via level set
    ///
    std::cout << "CholTriSolve:: solve matrix via TriSolveTools" << std::endl;

    timer.reset();        
#if defined(TACHO_USE_TRISOLVE_VARIANT)
    constexpr int variant = TACHO_USE_TRISOLVE_VARIANT;
#else
    constexpr int variant = 0;
#endif
    Tacho::TriSolveTools<value_type,scheduler_type,variant> 
      TS(t_perm, t_peri,
         N.getSupernodesInfo(),
         S.SupernodesTreeLevel(),
         nrhs);
    TS.initialize(device_level_cut, device_function_thres, verbose);
    TS.createStream(nstreams_solve);
    TS.prepareSolve(nstreams_prepare, verbose); 
    t = timer.seconds();    
    std::cout << "CholTriSolve:: TriSolve prepare::time = " << t << std::endl;

    timer.reset();        
    TS.solveCholesky(X, B, W, verbose); 
    t = timer.seconds();    
    std::cout << "CholTriSolve:: solve matrix::time = " << t << std::endl;
    TS.release(verbose);
#else
    ///
    /// solve via tasking
    ///
    std::cout << "CholTriSolve:: solve matrix" << std::endl;
    timer.reset();    
    if (serial) {
#if !defined (KOKKOS_ENABLE_CUDA)
      N.solveCholesky_Serial(X, B, W0, verbose);
#endif
    } else {
      N.solveCholesky_Parallel(X, B, W0, verbose);
    }
    t = timer.seconds();    
    std::cout << "CholTriSolve:: solve matrix::time = " << t << std::endl;
#endif

    const double res = N.computeRelativeResidual(X, B);
    std::cout << "CholTriSolve:: residual = " << res << std::endl;

    N.release(verbose);
  }
  Kokkos::finalize();

  return r_val;
}
