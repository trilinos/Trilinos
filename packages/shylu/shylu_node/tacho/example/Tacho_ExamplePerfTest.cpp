#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

#include "Tacho_CommandLineParser.hpp"

#ifdef TACHO_HAVE_VTUNE
#include "ittnotify.h"
#define TACHO_ITT_PAUSE __itt_pause()
#define TACHO_ITT_RESUME __itt_resume()
#else
#define TACHO_ITT_PAUSE 
#define TACHO_ITT_RESUME
#endif

#if defined( __INTEL_MKL__ )
#include "mkl_service.h"
#include "Tacho_Pardiso.hpp"
#endif

#if defined( TACHO_HAVE_SUITESPARSE )
#include "cholmod.h"
#endif

int main (int argc, char *argv[]) {
  int nthreads = 1; 

  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;
  int niter_solve = 50;

  int max_num_superblocks = 4;
  int mb = -1, nb = -1;
  int small_problem_thres = 1024;
  int front_update_mode = -1;

  bool test_tacho = true;
  bool test_pardiso = false;
  bool test_cholmod = false;

  bool use_same_ordering = true;

  Tacho::CommandLineParser opts("This is Tacho performance test comparing with Pardiso and Cholmod on OpenMP and Cuda spaces");

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
  opts.set_option<int>("front-update-mode", "Front update mode", &front_update_mode);

  // testing flags
  opts.set_option<bool>("test-tacho", "Flag for testing Tacho", &test_tacho);
#if defined( __INTEL_MKL__ )
  opts.set_option<bool>("test-pardiso", "Flag for testing Pardiso", &test_pardiso);
#endif
#if defined( TACHO_HAVE_SUITESPARSE )
  opts.set_option<bool>("test-cholmod", "Flag for testing Cholmod", &test_cholmod);
#endif

  opts.set_option<bool>("use-same-ordering", "Same Metis ordering is used for all tests", &use_same_ordering);

  TACHO_ITT_PAUSE;
  
  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; 

  Kokkos::initialize(argc, argv);
  Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);

  int r_val = 0;

  {
    /// basic typedef
    typedef int ordinal_type;
    typedef double value_type;

    /// crs matrix format and dense multi vector
    typedef Tacho::CrsMatrixBase<value_type,Kokkos::DefaultHostExecutionSpace> CrsMatrixBaseType;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,Kokkos::DefaultHostExecutionSpace> DenseMultiVectorType;
    //typedef Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace> OrdinalTypeArray;

    ///
    /// problem setting
    ///
    CrsMatrixBaseType A;
    {
      {
        std::ifstream in;
        in.open(file);
        if (!in.good()) {
          std::cout << "Failed in open the file: " << file << std::endl;
          return -1;
        }
      }
      Tacho::MatrixMarket<value_type>::read(file, A, verbose);
    }

    DenseMultiVectorType
      b("b", A.NumRows(), nrhs), // rhs multivector
      x("x", A.NumRows(), nrhs), // solution multivector
      t("t", A.NumRows(), nrhs); // temp workspace (store permuted rhs)

    Tacho::Graph graph(A.NumRows(), A.NumNonZeros(), A.RowPtr(), A.Cols());
#if   defined(TACHO_HAVE_METIS)
    Tacho::GraphTools_Metis G(graph);
#elif defined(TACHO_HAVE_SCOTCH)
    Tacho::GraphTools_Scotch G(graph);
#else
    Tacho::GraphTools_CAMD G(graph);
#endif
    G.reorder(verbose);

    {
      Tacho::Random<value_type> random;
      const ordinal_type m = A.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs)
        for (ordinal_type i=0;i<m;++i)
          b(i, rhs) = random.value();
    }

    // KNL - 34 tiles share 1MB L2 separate cache
    // Skylake - 33 MB shared L3
    // V100 - do I need to flush ?
    constexpr size_t LLC_CAPACITY = 34*1024*1024;
    Tacho::Flush<LLC_CAPACITY> flush;

    // -----------------------------------------------------------------
    if (test_pardiso) {
#if defined( __INTEL_MKL__ )
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
      if (use_same_ordering)
        pardiso.setParameter(5, 1); // user permutation is used for mkl permutation

      // somehow pardiso does not like symmetric full matrix (store only half)
      CrsMatrixBaseType Asym("Asym");
      Asym.createConfTo(A);
      {
        size_t nnz = 0;
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
                         (int*)G.PermVector().data(),
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
      {
        const auto peri = G.InvPermVector();
        const auto perm = G.PermVector();
        for (ordinal_type i=0;i<A.NumRows();++i) 
          peri(perm(i)) = i;  
      }

      TACHO_ITT_RESUME;
      r_val = pardiso.run(Pardiso::Factorize, 1);
      TACHO_ITT_PAUSE;

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
        printf("  Time (Multiple Solves) \n");                                                                                     
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

      const double res = Tacho::NumericTools<value_type,Kokkos::DefaultHostExecutionSpace>::computeRelativeResidual(A, x, b);
      std::cout << "PardisoChol:: residual = " << res << "\n\n";

      r_val = pardiso.run(Pardiso::ReleaseAll);
      if (r_val) {
        std::cout << "PardisoChol:: release error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::ReleaseAll) << std::endl;
      }
      std::cout << "\n\n";
#else
      std::cout << "MKL is NOT configured in Trilinos" << std::endl;
#endif
    }

    // -----------------------------------------------------------------
    if (test_cholmod) {
#if defined( TACHO_HAVE_SUITESPARSE )
      flush.run();

      Kokkos::Impl::Timer timer;
      double t_analyze = 0, t_factor = 0, t_solve = 0, t_solve_niter = 0;

      cholmod_sparse *AA ;
      cholmod_dense *xx, *bb, *rr;
      cholmod_factor *LL ;
      double one [2] = {1,0}, m1 [2] = {-1,0} ;    /* basic scalars */

      // cholmod handle
      cholmod_common c;
      cholmod_start (&c);    /* start CHOLMOD */
      c.nmethods = 1;
      c.method[0].ordering = use_same_ordering ? CHOLMOD_GIVEN : CHOLMOD_METIS;
      c.postorder = 1;

      AA = cholmod_allocate_sparse(A.NumRows(), A.NumRows(), A.NumNonZeros(),
                                   1, 1, 1, CHOLMOD_REAL, &c);
      {
        int *pp = (int*)AA->p;
        for (int i=0;i<(A.NumRows()+1);++i) 
          pp[i] = A.RowPtrBegin(i);
        
        int *ii = (int*)AA->i, jend = A.NumNonZeros();
        double *aa = (double*)AA->x;
        for (int j=0;j<jend;++j) {
          ii[j] = A.Col(j);
          aa[j] = A.Value(j);
        }
      }

      {
        const int m = AA->nrow, n = nrhs, lda = AA->nrow;
        bb = cholmod_allocate_dense(m, n, lda, CHOLMOD_REAL, &c);

        double *bbx = (double*)bb->x;
        for (int j=0;j<n;++j)
          for (int i=0;i<m;++i)
            bbx[i+j*lda] = b(i,j);
      }

      // reset memory usage to measure only for factors
      c.memory_inuse = 0;
      c.memory_usage = 0;

      timer.reset();
      if (use_same_ordering)
        LL = cholmod_analyze_p (AA, G.PermVector().data(), NULL, 0, &c) ;    /* analyze */
      else
        LL = cholmod_analyze (AA, &c) ;    /* analyze */
      t_analyze = timer.seconds();

      if (verbose) {
        printf("CHOLMOD: Analyze\n");
        printf("================\n");
        printf("  Time\n");
        printf("             total time spent:                                 %10.6f s\n", t_analyze);
        printf("\n");
        printf("  Linear system A\n");
        printf("             number of equations:                              %10zu\n", AA->nrow);
        printf("             max number of nonzeros:                           %10zu\n", AA->nzmax);
        printf("\n");        
        printf("  Factors L\n");
        printf("             number of nonzeros:                               %10.0f\n", c.lnz);
        //printf("             current method:                                   %10d\n", c.current);
        printf("\n");                
      }

      TACHO_ITT_RESUME;
      timer.reset();
      cholmod_factorize (AA, LL, &c) ;    /* factorize */
      t_factor = timer.seconds();
      TACHO_ITT_PAUSE;

      if (verbose) {
        printf("CHOLMOD: Factorize\n");
        printf("==================\n");
        printf("  Time\n");
        printf("             total time spent:                                 %10.6f s\n", t_factor);
        printf("\n");                        
        printf("  Property\n");
        printf("             is_super:                                         %10d\n", LL->is_super);
        printf("             ordering:                                         %10d\n", LL->ordering);
        printf("               0 - natural\n");
        printf("               1 - given\n");
        printf("               2 - AMD\n");
        printf("               3 - METIS\n");
        printf("               4 - NESDIS (CHOLMOD nd)\n");
        printf("               5 - AMD for A, COLAMD for A*A'\n");
        printf("\n");                
        printf("  Memory\n");
        printf("             memory used:                                      %10.6f MB\n", double(c.memory_inuse)/1024/1024);
        printf("             peak memory used in factorization:                %10.6f MB\n", double(c.memory_usage)/1024/1024);
        printf("\n");                
        printf("  FLOPs\n");
        printf("             gflop   for numeric factorization:                %10.6f GFLOP\n", c.fl/1024/1024/1024);
        printf("             gflop/s for numeric factorization:                %10.6f GFLOP/s\n", c.fl/1024/1024/1024/t_factor);
        printf("\n");               
      }

      timer.reset();
      for (int iter=0;iter<niter_solve;++iter) {
        xx = cholmod_solve (CHOLMOD_A, LL, bb, &c) ;    /* solve Ax=b */
      }
      t_solve_niter = timer.seconds();
      t_solve = t_solve_niter / double(niter_solve);

      if (verbose) {
        printf("CHOLMOD: Solve\n");
        printf("==============\n");
        printf("  Time (Multiple Solves) \n");
        printf("             total time spent for %3d numeric solve:           %10.6f s\n", niter_solve, t_solve_niter);
        printf("             average time spent for a single numeric solve:    %10.6f s\n", t_solve);
        printf("\n\n");
      }

      if (test_tacho) {
        double diff = 0, norm = 0;
        const int m = AA->nrow, n = nrhs, lda = AA->nrow;
        double *xxx = (double*)xx->x;
        for (int j=0;j<n;++j) 
          for (int i=0;i<m;++i) {
            norm += x(i,j)*x(i,j);
            const double tmp = xxx[i+j*lda] - x(i,j);
            diff += tmp*tmp;
          }
        printf ("CHOLMOD: diff to tacho %10.6e\n", diff/norm);
      }

      rr = cholmod_copy_dense (bb, &c) ;    /* r = b */
      cholmod_sdmult (AA, 0, m1, one, xx, rr, &c) ;    /* r = r-Ax */
      printf ("CHOLMOD: residual %10.6e\n", cholmod_norm_dense (rr, 0, &c));

      cholmod_free_factor (&LL, &c) ;    /* free matrices */
      cholmod_free_sparse (&AA, &c) ;
      cholmod_free_dense (&rr, &c) ;
      cholmod_free_dense (&xx, &c) ;
      cholmod_free_dense (&bb, &c) ;
      cholmod_finish (&c) ; /* finish CHOLMOD */
      printf("CHOLMOD: Finished\n");
      printf("=================\n");
#else
      std::cout << "CHOLMOD is NOT configured in Trilinos" << std::endl;
#endif
    }

    // -----------------------------------------------------------------
    if (test_tacho) {
      flush.run();

      Kokkos::Impl::Timer timer;
      double t_solve = 0, t_solve_niter = 0;

      ///
      /// tacho
      ///
      Tacho::Solver<value_type,Kokkos::DefaultHostExecutionSpace> solver;

      //solver.setMatrixType(sym, posdef);
      solver.setVerbose(verbose);
      solver.setMaxNumberOfSuperblocks(max_num_superblocks);
      solver.setSmallProblemThresholdsize(small_problem_thres);
      solver.setBlocksize(mb);
      solver.setPanelsize(nb);
      solver.setFrontUpdateMode(front_update_mode);

      /// inputs are used for graph reordering and analysis
      solver.analyze(A.NumRows(),
                     A.RowPtr(),
                     A.Cols(),
                     G.PermVector(),
                     G.InvPermVector());

      /// symbolic structure can be reused
      TACHO_ITT_RESUME;
      solver.factorize(A.Values());
      TACHO_ITT_PAUSE;

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


  }
  Kokkos::finalize();
  return r_val;
}

