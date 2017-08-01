#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_Graph.hpp"
#include "TachoExp_SymbolicTools.hpp"

#if defined(HAVE_SHYLUTACHO_SCOTCH)
#include "TachoExp_GraphTools_Scotch.hpp"
#endif

#if defined(HAVE_SHYLUTACHO_METIS)
#include "TachoExp_GraphTools_Metis.hpp"
#endif

#include "TachoExp_GraphTools_CAMD.hpp"

#include "TachoExp_NumericTools.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif

using namespace Tacho;
using namespace Tacho::Experimental;

int main (int argc, char *argv[]) {
  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Pardiso Chol algorithms on Kokkos::OpenMP execution space.\n");

  bool serial = false;
  clp.setOption("enable-serial", "disable-verbose", &serial, "Flag for invoking serial algorithm");

  int nthreads = 1;
  clp.setOption("kokkos-threads", &nthreads, "Number of threads");

  bool verbose = true;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  std::string file = "test.mtx";
  clp.setOption("file", &file, "Input file (MatrixMarket SPD matrix)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of RHS vectors");

  int serial_thres_size = -1; // 32 is better
  clp.setOption("serial-thres", &serial_thres_size, "Serial threshold");  
  
  int mb = 0;
  clp.setOption("mb", &mb, "Blocksize");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  //if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

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
#if   defined(HAVE_SHYLUTACHO_METIS)
    GraphTools_Metis T(G);
#elif defined(HAVE_SHYLUTACHO_SCOTCH)
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
