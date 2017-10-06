#include "ShyLU_NodeTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_Graph.hpp"
#include "TachoExp_SymbolicTools.hpp"

#if defined(TACHO_HAVE_SCOTCH)
#include "TachoExp_GraphTools_Scotch.hpp"
#endif

#if defined(TACHO_HAVE_METIS)
#include "TachoExp_GraphTools_Metis.hpp"
#endif

#include "TachoExp_GraphTools_CAMD.hpp"

#include "TachoExp_NumericTools.hpp"

#include "TachoExp_CommandLineParser.hpp"

#ifdef TACHO_HAVE_MKL
#include "mkl_service.h"
#endif

using namespace Tacho;
using namespace Tacho::Experimental;

int main (int argc, char *argv[]) {
  CommandLineParser opts("This example program measure the performance of Tacho on Kokkos::OpenMP");

  bool serial = false;
  int nthreads = 1;
  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;
  int serial_thres_size = -1; // 32 is better
  int mb = 0;

  opts.set_option<bool>("serial", "Flag to use serial algorithm", &serial);
  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("serial-thres", "Serialization threshold size", &serial_thres_size);
  opts.set_option<int>("mb", "Blocksize", &mb);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);
  if (std::is_same<Kokkos::DefaultHostExecutionSpace,Kokkos::Serial>::value) 
    std::cout << "Kokkos::Serial\n";
  else
    Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);
  
  int r_val = 0;
  
  {
    typedef double value_type;
    typedef CrsMatrixBase<value_type> CrsMatrixBaseType;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,Kokkos::DefaultHostExecutionSpace> DenseMatrixBaseType;
    
    Kokkos::Impl::Timer timer;
    double t = 0.0;

    std::cout << "CholSupernodes:: import input file = " << file << std::endl;
    CrsMatrixBaseType A("A");
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
      A = MatrixMarket<value_type>::read(file);
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

    NumericTools<value_type,Kokkos::DefaultHostExecutionSpace> 
      N(A.NumRows(), A.RowPtr(), A.Cols(), // A.Values(),
        T.PermVector(), T.InvPermVector(),
        S.NumSupernodes(), S.Supernodes(),
        S.gidSuperPanelPtr(), S.gidSuperPanelColIdx(),
        S.sidSuperPanelPtr(), S.sidSuperPanelColIdx(), S.blkSuperPanelColIdx(),
        S.SupernodesTreeParent(), S.SupernodesTreePtr(), S.SupernodesTreeChildren(), S.SupernodesTreeRoots());
    N.setSerialThresholdSize(serial_thres_size);

    std::cout << "CholSupernodes:: factorize matrix" << std::endl;
    timer.reset();    
    if (serial) {
      N.factorizeCholesky_Serial(A.Values(), verbose);
    } else {
      if (mb > 0) 
        N.factorizeCholesky_ParallelByBlocks(A.Values(), mb, verbose);
      else
        N.factorizeCholesky_Parallel(A.Values(), verbose);
    }
    t = timer.seconds();    
    std::cout << "CholSupernodes:: factorize matrix::time = " << t << std::endl;
    
    DenseMatrixBaseType 
      B("B", A.NumRows(), nrhs), 
      X("X", A.NumRows(), nrhs), 
      Y("Y", A.NumRows(), nrhs);
    
    {
      Random<value_type> random;
      const ordinal_type m = A.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs)
        for (ordinal_type i=0;i<m;++i) 
          B(i, rhs) = random.value();
    }

    std::cout << "CholSupernodes:: solve matrix" << std::endl;
    timer.reset();    
    if (serial) {
      N.solveCholesky_Serial(X, B, Y, verbose);
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
