#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

#include "TachoExp_CommandLineParser.hpp"

#if defined( TACHO_HAVE_MKL )
#include "mkl_service.h"
#include "Tacho_Pardiso.hpp"
#endif

#if defined( TACHO_HAVE_SUITESPARSE )
#endif

int main (int argc, char *argv[]) {
  int nthreads = 1; 

  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;
  int niter_solve = 50;

  int max_num_superblocks = 4;
  int mb = -1, nb = -1;

  bool test_tacho = true;
  bool test_pardiso = false;
  bool test_cholmod = false;

  Tacho::Experimental::CommandLineParser opts("This is Tacho performance test comparing with Pardiso and Cholmod on OpenMP and Cuda spaces");

  // threading environment
  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);

  // common testing environment
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("niter-solve", "Number of solve iterations for timing", &niter_solve);

  // tacho parameters
  opts.set_option<int>("max-num-superblocks", "Max number of superblocks", &max_num_superblocks);
  opts.set_option<int>("mb", "Internal block size", &mb);
  opts.set_option<int>("nb", "Internal panel size", &nb);

  // testing flags
  opts.set_option<bool>("test-tacho", "Flag for testing Tacho", &test_tacho);
  opts.set_option<bool>("test-pardiso", "Flag for testing Pardiso", &test_pardiso);
  opts.set_option<bool>("test-cholmod", "Flag for testing Cholmod", &test_cholmod);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; 

  Kokkos::initialize(argc, argv);
  if (std::is_same<Kokkos::DefaultHostExecutionSpace,Kokkos::Serial>::value)
    std::cout << "Kokkos::Serial\n";
  else
    Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);

  int r_val = 0;

  {
    /// basic typedef
    typedef int ordinal_type;
    typedef size_t size_type;
    typedef double value_type;

    /// crs matrix format and dense multi vector
    typedef Tacho::Experimental::CrsMatrixBase<value_type,Kokkos::DefaultHostExecutionSpace> CrsMatrixBaseType;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,Kokkos::DefaultHostExecutionSpace> DenseMultiVectorType;
    typedef Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace> OrdinalTypeArray;

    ///
    /// problem setting
    ///
    CrsMatrixBaseType A("A");
    {
      {
        std::ifstream in;
        in.open(file);
        if (!in.good()) {
          std::cout << "Failed in open the file: " << file << std::endl;
          return -1;
        }
      }
      A = Tacho::Experimental::MatrixMarket<value_type>::read(file, verbose);
    }

    DenseMultiVectorType
      b("b", A.NumRows(), nrhs), // rhs multivector
      x("x", A.NumRows(), nrhs), // solution multivector
      t("t", A.NumRows(), nrhs); // temp workspace (store permuted rhs)

    OrdinalTypeArray perm("perm", A.NumRows()), peri("peri", A.NumRows());

    {
      Tacho::Experimental::Random<value_type> random;
      const ordinal_type m = A.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs)
        for (ordinal_type i=0;i<m;++i)
          b(i, rhs) = random.value();
    }

    // 34 tiles share 1MB L2 separate cache
    constexpr size_t LLC_CAPACITY = 34*1024*1024;
    Tacho::Experimental::Flush<LLC_CAPACITY> flush;


    // -----------------------------------------------------------------
    if (test_tacho) {
      flush.run();

      Kokkos::Impl::Timer timer;
      double t_solve = 0, t_solve_niter = 0;

      ///
      /// tacho
      ///
      Tacho::Solver<value_type,Kokkos::DefaultHostExecutionSpace> solver;
      solver.setVerbose(verbose);
      solver.setMaxNumberOfSuperblocks(max_num_superblocks);

      /// inputs are used for graph reordering and analysis
      solver.analyze(A.NumRows(),
                     A.RowPtr(),
                     A.Cols());

      /// symbolic structure can be reused
      solver.factorize(A.Values());

      solver.setVerbose(0); // disable verbose out for the iteration
      timer.reset();
      for (int iter=0;iter<niter_solve;++iter) {
        solver.solve(x, b, t);
      }
      t_solve_niter = timer.seconds();
      t_solve = t_solve_niter / double(niter_solve);

      solver.setVerbose(verbose);
      solver.solve(x, b, t);

      if (verbose) {
        printf("  Time (Multiple Solves) \n");                                                                                     
        printf("             total time spent for %3d numeric solve:           %10.6f s\n", niter_solve, t_solve_niter);
        printf("             average time spent for a single numeric solve:    %10.6f s\n", t_solve);
        printf("\n\n");
      }
      
      const double res = solver.computeRelativeResidual(A.Values(), x, b);
      std::cout << "TachoSolver: residual = " << res << "\n\n";
      
      solver.release();
    }

    // -----------------------------------------------------------------
    if (test_pardiso) {
#if defined( TACHO_HAVE_MKL )
      flush.run();

      Kokkos::Impl::Timer timer;
      double t_solve = 0, t_solve_niter = 0;

      typedef Tacho::Pardiso Pardiso;

      // mkl nthreads setting
      // mkl_set_dynamic(0);
      // mkl_set_num_threads(nthreads);

      Pardiso pardiso;

      constexpr int AlgoChol = 2;
      r_val = pardiso.init<value_type,AlgoChol>();
      if (r_val) {
        std::cout << "PardisoChol:: Pardiso init error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      }
      pardiso.setParameter(2, 2); // metis ordering is used 
      //pardiso.setParameter(5, 1); // user permutation is used for mkl permutation
      //pardiso.setParameter(5, 2); // mkl permutation is written on perm

      // somehow pardiso does not like symmetric full matrix (store only half)
      CrsMatrixBaseType Asym("Asym");
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

      // 32bit vs 64bit integers; A uses size_t for size array
      Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace> rowptr("rowptr", Asym.NumRows()+1);
      for (ordinal_type i=0;i<=Asym.NumRows();++i)
        rowptr(i) = Asym.RowPtrBegin(i);

      pardiso.setProblem(Asym.NumRows(),
                         (double*)Asym.Values().data(),
                         (int*)rowptr.data(),// (int*)Asym.RowPtr().data(),
                         (int*)Asym.Cols().data(),
                         (int*)perm.data(),
                         nrhs,
                         (double*)b.data(),
                         (double*)x.data());
      
      r_val = pardiso.run(Pardiso::Analyze, 1);
      if (r_val) {                                                    
        std::cout << "PardisoChol:: Pardiso analyze error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;                 
      } else {                                                         
        pardiso.showStat(std::cout, Pardiso::Analyze) << std::endl;   
      }

      // compute inverse permutation
      for (ordinal_type i=0;i<A.NumRows();++i) 
        peri(perm(i)) = i;  

      r_val = pardiso.run(Pardiso::Factorize, 1);
      if (r_val) {                                                    
        std::cout << "PardisoChol:: Pardiso factorize error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;                 
      } else {                                                         
        pardiso.showStat(std::cout, Pardiso::Factorize) << std::endl;   
      }

      timer.reset();
      for (int iter=0;iter<niter_solve;++iter) {
        r_val = pardiso.run(Pardiso::Solve, 0);
      }
      t_solve_niter = timer.seconds();
      t_solve = t_solve_niter / double(niter_solve);

      r_val = pardiso.run(Pardiso::Solve, 1);

      if (verbose) {
        printf("  Time\n");                                                                                     
        printf("             total time spent for %3d numeric solve:           %10.6f s\n", niter_solve, t_solve_niter);
        printf("             average time spent for a single numeric solve:    %10.6f s\n", t_solve);
        printf("\n\n");
      }

      if (r_val) {                                                    
        std::cout << "PardisoChol:: Pardiso solve error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;                 
      } else {                                                         
        pardiso.showStat(std::cout, Pardiso::Solve) << std::endl;   
      }

      const double res = Tacho::Experimental::NumericTools<value_type,Kokkos::DefaultHostExecutionSpace>::computeRelativeResidual(A, x, b);
      std::cout << "PardisoChol:: residual = " << res << "\n\n";

      r_val = pardiso.run(Pardiso::ReleaseAll);
      if (r_val) {
        std::cout << "PardisoChol:: release error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::ReleaseAll) << std::endl;
      }
#else
      cout << "MKL is NOT configured in Trilinos" << endl;
#endif
    }
  }
  Kokkos::finalize();
  return r_val;
}


