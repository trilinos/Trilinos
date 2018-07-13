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

#include "Tacho_CommandLineParser.hpp"

#ifdef TACHO_HAVE_MKL
#include "mkl_service.h"
#endif

using namespace Tacho;

int main (int argc, char *argv[]) {
  CommandLineParser opts("This example program measure the performance of Tacho on Kokkos::OpenMP");

  bool serial = false;
  int nthreads = 1;
  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;
  int serial_thres_size = -1; // 32 is better
  int mb = 0;
  int nb = 0;

  opts.set_option<bool>("serial", "Flag to use serial algorithm", &serial);
  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("serial-thres", "Serialization threshold size", &serial_thres_size);
  opts.set_option<int>("mb", "Blocksize", &mb);
  opts.set_option<int>("nb", "Panelsize", &nb);

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

  typedef Kokkos::DefaultExecutionSpace exec_space;
  //typedef Kokkos::DefaultHostExecutionSpace exec_space;
  typedef Kokkos::DefaultHostExecutionSpace host_space;

  printExecSpaceConfiguration<exec_space> ("DeviceSpace", false);
  printExecSpaceConfiguration<host_space> ("HostSpace",   false);
  
  int r_val = 0;
  
  {
    typedef double value_type;
    typedef CrsMatrixBase<value_type,host_space> CrsMatrixBaseTypeHost;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,exec_space> DenseMatrixBaseType;
    
    Kokkos::Impl::Timer timer;
    double t = 0.0;

    std::cout << "CholSupernodes:: import input file = " << file << std::endl;
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
      MatrixMarket<value_type>::read(file, A);
    }
    Graph G(A);
    t = timer.seconds();
    std::cout << "CholSupernodes:: import input file::time = " << t << std::endl;

    std::cout << "CholSupernodes:: analyze matrix" << std::endl;
    timer.reset();
#if   defined(TACHO_HAVE_METIS)
    GraphTools_Metis T(G);
#elif defined(TACHO_HAVE_SCOTCH)
    GraphTools_Scotch T(G);
#else
    GraphTools_CAMD T(G);
#endif
    T.reorder(verbose);
    
    SymbolicTools S(A, T);
    S.symbolicFactorize(verbose);
    t = timer.seconds();
    std::cout << "CholSupernodes:: analyze matrix::time = " << t << std::endl;

    typedef typename exec_space::memory_space device_memory_space; 
    
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
    
    NumericTools<value_type,exec_space> 
      N(A.NumRows(), a_row_ptr, a_cols,
        t_perm, t_peri,
        S.NumSupernodes(), s_supernodes,
        s_gid_spanel_ptr, s_gid_spanel_colidx,
        s_sid_spanel_ptr, s_sid_spanel_colidx, s_blk_spanel_colidx,
        s_snodes_tree_parent, s_snodes_tree_ptr, s_snodes_tree_children,
        S.SupernodesTreeRoots());

    N.setSerialThresholdSize(serial_thres_size);

    std::cout << "CholSupernodes:: factorize matrix" << std::endl;
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
    std::cout << "CholSupernodes:: factorize matrix::time = " << t << std::endl;
    
    DenseMatrixBaseType 
      B("B", A.NumRows(), nrhs), 
      X("X", A.NumRows(), nrhs), 
      Y("Y", A.NumRows(), nrhs);

    {
      Kokkos::Random_XorShift64_Pool<exec_space> random(13718);      
      Kokkos::fill_random(B, random, value_type(1));
    }

    std::cout << "CholSupernodes:: solve matrix" << std::endl;
    timer.reset();    
    if (serial) {
#if !defined (KOKKOS_ENABLE_CUDA)
      N.solveCholesky_Serial(X, B, Y, verbose);
#endif
    } else {
      N.solveCholesky_Parallel(X, B, Y, verbose);
    }
    t = timer.seconds();    
    std::cout << "CholSupernodes:: solve matrix::time = " << t << std::endl;

    const double res = N.computeRelativeResidual(X, B);
    //const double eps = std::numeric_limits<double>::epsilon()*100;
    std::cout << "CholSupernodes:: residual = " << res << std::endl;
  }
  Kokkos::finalize();

  return r_val;
}
