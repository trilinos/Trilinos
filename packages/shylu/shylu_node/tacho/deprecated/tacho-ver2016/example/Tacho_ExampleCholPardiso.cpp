#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLU_NodeTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_MatrixMarket.hpp"

#ifdef HAVE_SHYLU_NODETACHO_MKL
#include "mkl_service.h"
#include "Tacho_ExamplePardiso.hpp"
#endif

using namespace Tacho;
using namespace Tacho::Experimental;

int main (int argc, char *argv[]) {
  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Pardiso Chol algorithms on Kokkos::OpenMP execution space.\n");

  int nthreads = 1;
  clp.setOption("kokkos-threads", &nthreads, "Number of threads");

  bool verbose = false;
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
#ifdef HAVE_SHYLU_NODETACHO_MKL
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
    CrsMatrixBaseType A("A");
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
      
      // somehow pardiso does not like symmetric matrix
      CrsMatrixBaseType Atmp("Atmp");
      Atmp.createConfTo(A);
      {
        size_type nnz = 0;
        for (ordinal_type i=0;i<A.NumRows();++i) {
          Atmp.RowPtrBegin(i) = nnz;
          for (ordinal_type idx=A.RowPtrBegin(i);idx<A.RowPtrEnd(i);++idx) {
            if (i <= A.Col(idx)) {
              Atmp.Col(nnz) = A.Col(idx);
              Atmp.Value(nnz) = A.Value(idx);
              ++nnz;
            }
          }
          Atmp.RowPtrEnd(i) = nnz;
        }
      }
      A = Atmp;
    }
    t = timer.seconds();

    // 32bit vs 64bit integers; A uses size_t for size array
    Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace> rowptr("rowptr", A.NumRows()+1);
    {      
      for (ordinal_type i=0;i<=A.NumRows();++i)
        rowptr(i) = A.RowPtrBegin(i);
    }    
    std::cout << "PardisoChol:: import input file::time = " << t << std::endl;
    
    DenseMatrixBaseType 
      BB("BB", A.NumRows(), nrhs), 
      XX("XX", A.NumRows(), nrhs), 
      RR("RR", A.NumRows(), nrhs),
      PP("PP",  A.NumRows(), 1);
    
    {
      const auto m = A.NumRows();
      srand(time(NULL));
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) 
          XX(i, rhs) = ((value_type)rand()/(RAND_MAX));
        
        // matvec
        Kokkos::DefaultHostExecutionSpace::fence();
        Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, m),
                             [&](const ordinal_type i) {
                               value_type tmp = 0;
                               for (ordinal_type j=A.RowPtrBegin(i);j<A.RowPtrEnd(i);++j)
                                 tmp += A.Value(j)*XX(A.Col(j), rhs);
                               BB(i, rhs) = tmp;
                             } );
        Kokkos::DefaultHostExecutionSpace::fence();
      }
      Kokkos::deep_copy(RR, XX);
    }

    pardiso.setProblem(A.NumRows(),
                       (double*)A.Values().data(),
                       (int*)rowptr.data(),// (int*)A.RowPtr().data(),
                       (int*)A.Cols().data(),
                       (int*)PP.data(),
                       nrhs,
                       (double*)BB.data(),
                       (double*)XX.data());
    
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
      double error = 0, norm = 0;
      const auto m = A.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) {
          {
            const auto val = std::abs(XX(i, rhs) - RR(i, rhs));
            error += val*val;
          }
          {
            const auto val = std::abs(RR(i, rhs));
            norm  += val*val;
          }
        }
      }
      std::cout << "PardisoChol:: error = " << error << " , norm = " << norm << std::endl;
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
