#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#include "TachoExp_NumericTools.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#include "Tacho_Pardiso.hpp"
#endif

using namespace Tacho;
using namespace Tacho::Experimental;

int main (int argc, char *argv[]) {
  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Pardiso Chol algorithms on Kokkos::OpenMP execution space.\n");

  int nthreads = 1;
  clp.setOption("kokkos-threads", &nthreads, "Number of threads");

  bool verbose = true;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  std::string file_input = "test.mtx";
  clp.setOption("file", &file_input, "Input file (MatrixMarket SPD matrix)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of RHS vectors");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  //if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  const bool skip_factorize = false, skip_solve = false;
  
  Kokkos::initialize(argc, argv);
  Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);

  int r_val = 0;
#ifdef HAVE_SHYLUTACHO_MKL
  {
    typedef double value_type;
    typedef CrsMatrixBase<value_type> CrsMatrixBaseType;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,Kokkos::DefaultHostExecutionSpace> DenseMatrixBaseType;
    
    // mkl nthreads setting 
    mkl_set_dynamic(0);
    mkl_set_num_threads(nthreads);
    
    Kokkos::Impl::Timer timer;
    double t = 0.0;
    Pardiso pardiso;

    constexpr int AlgoChol = 2;
    std::cout << "PardisoChol:: init" << std::endl;
    {
      timer.reset();
      r_val = pardiso.init<value_type,AlgoChol>();
      t = timer.seconds();
      
      if (r_val) {
        std::cout << "PardisoChol:: Pardiso init error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      }
    }
    std::cout << "PardisoChol:: init ::time = " << t << std::endl;
    
    std::cout << "PardisoChol:: import input file = " << file_input << std::endl;
    CrsMatrixBaseType A("A"), Asym("Asym");
    timer.reset();
    {
      {
        std::ifstream in;
        in.open(file_input);
        if (!in.good()) {
          std::cout << "Failed in open the file: " << file_input << std::endl;
          return -1;
        }
      }
      A = MatrixMarket<value_type>::read(file_input);
      
      // somehow pardiso does not like symmetric full matrix (store only half)
      Asym.createConfTo(A);
      {
        size_type nnz = 0;
        for (ordinal_type i=0;i<A.NumRows();++i) {
          Asym.RowPtrBegin(i) = nnz;
          for (ordinal_type idx=A.RowPtrBegin(i);idx<A.RowPtrEnd(i);++idx) {
            if (i <= A.Col(idx)) {
              Asym.Col(nnz) = A.Col(idx);
              Asym.Value(nnz) = A.Value(idx);
              ++nnz;
            }
          }
          Asym.RowPtrEnd(i) = nnz;
        }
      }
    }
    t = timer.seconds();

    // 32bit vs 64bit integers; A uses size_t for size array
    Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace> rowptr("rowptr", Asym.NumRows()+1);
    {      
      for (ordinal_type i=0;i<=Asym.NumRows();++i)
        rowptr(i) = Asym.RowPtrBegin(i);
    }    
    std::cout << "PardisoChol:: import input file::time = " << t << std::endl;
    
    DenseMatrixBaseType 
      B("B", Asym.NumRows(), nrhs), 
      X("X", Asym.NumRows(), nrhs), 
      P("P", Asym.NumRows(), 1);
    
    {
      const auto m = Asym.NumRows();
      Random<value_type> random;
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) 
        for (ordinal_type i=0;i<m;++i) 
          B(i, rhs) = random.value();
      Kokkos::deep_copy(X, B);
    }
    
    pardiso.setProblem(Asym.NumRows(),
                       (double*)Asym.Values().data(),
                       (int*)rowptr.data(),// (int*)Asym.RowPtr().data(),
                       (int*)Asym.Cols().data(),
                       (int*)P.data(),
                       nrhs,
                       (double*)B.data(),
                       (double*)X.data());
    
    std::cout << "PardisoChol:: analyze matrix" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::Analyze);
      t = timer.seconds();
      
      if (r_val) {
        std::cout << "PardisoChol:: Pardiso analyze error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::Analyze) << std::endl;
      }
    }
    std::cout << "PardisoChol:: analyze matrix::time = " << t << std::endl;

    if (!skip_factorize) {
      std::cout << "PardisoChol:: factorize matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Factorize);
        t = timer.seconds();
        
        if (r_val) {
          std::cout << "PardisoChol:: Pardiso factorize error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Factorize) << std::endl;
        }
      }
      std::cout << "PardisoChol:: factorize matrix::time = " << t << std::endl;
    }

    if (!skip_factorize && !skip_solve) {
      std::cout << "PardisoChol:: solve matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Solve);
        t = timer.seconds();
        
        if (r_val) {
          std::cout << "PardisoChol:: Pardiso solve error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Solve) << std::endl;
        }
      }
      std::cout << "PardisoChol:: solve matrix::time = " << t << std::endl;
    }
    
    {
      const double res = NumericTools<value_type,Kokkos::DefaultHostExecutionSpace>::computeRelativeResidual(A, X, B);
      std::cout << "PardisoChol:: residual = " << res << std::endl;
    }
    
    std::cout << "PardisoChol:: release all" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::ReleaseAll);
      t = timer.seconds();
      
      if (r_val) {
        std::cout << "PardisoChol:: release error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::ReleaseAll) << std::endl;
      }
    }
    std::cout << "PardisoChol:: release all::time = " << t << std::endl;
  }
#else
  r_val = -1;
  cout << "MKL is NOT configured in Trilinos" << endl;
#endif
  
  Kokkos::finalize();

  return r_val;
}
