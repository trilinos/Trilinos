#include "ShyLU_NodeTacho_config.h"

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho.hpp"
#include "Tacho_CommandLineParser.hpp"
#include "Tacho_Solver.hpp"

#if defined (KOKKOS_ENABLE_CUDA)
#include "Tacho_CuSparseTriSolve.hpp"
#endif

using namespace Tacho;

int main (int argc, char *argv[]) {

  CommandLineParser opts("This example program measure the performance of cuSparse TriSolve on Kokkos::Cuda");

  /// tacho factorization
  int max_num_superblocks = 32;
  int sym = 3;
  int posdef = 1;
  int small_problem_thres = 1024;
  int mb = 64;
  int nb = 64;

  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;

  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);

  opts.set_option<int>("max-num-superblocks", "Max number of superblocks", &max_num_superblocks);
  opts.set_option<int>("symmetric", "Symmetric type: 0 - unsym, 1 - structure sym, 2 - symmetric, 3 - hermitian", &sym);
  opts.set_option<int>("posdef", "Positive definite: 0 - indef, 1 - positive definite", &posdef);
  opts.set_option<int>("small-problem-thres", "LAPACK is used smaller than this thres", &small_problem_thres);
  opts.set_option<int>("mb", "Internal block size", &mb);
  opts.set_option<int>("nb", "Internal panel size", &nb);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);

  const bool detail = false;

  typedef double value_type;

  typedef UseThisDevice<Kokkos::Cuda>::device_type device_type;
  typedef UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;
  
  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace",   detail);

  typedef Kokkos::TaskSchedulerMultiple<typename device_type::execution_space> scheduler_type;

  Kokkos::Impl::Timer timer;
  int r_val = 0;
#if defined(KOKKOS_ENABLE_CUDA)
  {
    ///
    /// read from crs matrix
    ///
    typedef Tacho::CrsMatrixBase<value_type,host_device_type> CrsMatrixBaseTypeHost;    
    typedef Tacho::CrsMatrixBase<value_type,device_type> CrsMatrixBaseType;    
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> DenseMultiVectorType;

    /// read a spd matrix of matrix market format
    CrsMatrixBaseTypeHost h_A;
    {
      std::ifstream in;
      in.open(file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file << std::endl;
        return -1;
      }
      Tacho::MatrixMarket<value_type>::read(file, h_A, verbose);
    }
    
    ///
    /// A on device
    ///
    CrsMatrixBaseType A;
    A.createConfTo(h_A);
    A.copy(h_A);
    
    ///
    /// Tacho Solver (factorization)
    ///
    Tacho::Solver<value_type,scheduler_type> solver;
    {
      solver.setMatrixType(sym, posdef);
      solver.setVerbose(verbose);
      solver.setMaxNumberOfSuperblocks(max_num_superblocks);
      solver.setSmallProblemThresholdsize(small_problem_thres);
      solver.setBlocksize(mb);
      solver.setPanelsize(nb);
    }

    /// inputs are used for graph reordering and analysis
    solver.analyze(A.NumRows(),
                   A.RowPtr(),
                   A.Cols());
    
    /// symbolic structure can be reused
    solver.factorize(A.Values());

    /// crs matrix form of factors
    CrsMatrixBaseType F;
    {
      timer.reset();
      solver.exportFactorsToCrsMatrix(F);
      const double t_factor_export = timer.seconds();
      if (verbose) {
        printf("ExampleCuSparseTriSolve: Construction of CrsMatrix of factors\n");
        printf("=============================================================\n");
        printf("  Time\n");
        printf("             time for construction of F in CRS:               %10.6f s\n", t_factor_export); 
        printf("\n");

        F.showMe(std::cout, false);
      }      
    }

    ///
    /// CuSparseTriSolve
    ///
    CuSparseTriSolve trisolve;
    trisolve.setVerbose(verbose);
    
    ///
    /// CuSparse analyze
    ///
    {
      trisolve.analyze(F.NumRows(), 
                       F.RowPtr(),
                       F.Cols(),
                       F.Values());
    }

    ///
    /// random right hand side
    ///
    DenseMultiVectorType 
      b("b", F.NumRows(), nrhs), // rhs multivector
      x("x", F.NumRows(), nrhs), // solution multivector
      t("t", F.NumRows(), nrhs), // temporary workvector
      bb("bb", F.NumRows(), nrhs), // temp workspace (store permuted rhs)
      xx("xx", F.NumRows(), nrhs); // temp workspace (store permuted rhs)
    
    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);
      Kokkos::fill_random(b, random, value_type(1));
    }

    ///
    /// solve
    ///
    {
      const auto perm = solver.getPermutationVector();
      const auto peri = solver.getInversePermutationVector();

      const int niter = 1;
#if 0
      timer.reset();
      
      ///
      /// for now, using cuda graph does not gain performance
      /// for the problem size of our interest
      ///


      /// later we can create a permutation to be included 
      /// in the cuda graph
      trisolve.solve_capture(xx, bb, t);
      const double t_solve_capture = timer.seconds();

      timer.reset();
      for (int iter=0;iter<niter;++iter) {
        /// in order to use cuda graph RHS buffers should not be 
        /// reallocated
        applyRowPermutationToDenseMatrix(bb, b, perm);      
        trisolve.solve_launch();
        applyRowPermutationToDenseMatrix(x, xx, peri);            
        Kokkos::fence();
      }
      const double t_solve = timer.seconds();
      if (verbose) {
        printf("ExampleCuSparseTriSolve: P b, solve, and P^{-1} x\n");
        printf("=================================================\n");
        printf("  Time\n");
        printf("             time for capture:                                %10.6f s\n", t_solve_capture); 
        printf("             time for permute and solve:                      %10.6f s\n", t_solve); 
        printf("\n");
      }
#else  
      timer.reset();
      for (int iter=0;iter<niter;++iter) {
        applyRowPermutationToDenseMatrix(bb, b, perm);      
        trisolve.solve(xx, bb, t);
        applyRowPermutationToDenseMatrix(x, xx, peri);            
        Kokkos::fence();
      }
      const double t_solve = timer.seconds();
      if (verbose) {
        printf("ExampleCuSparseTriSolve: P b, solve, and P^{-1} x\n");
        printf("=================================================\n");
        printf("  Time\n");
        printf("             time for permute and solve:                      %10.6f s\n", t_solve); 
        printf("\n");
      }
#endif
    }

    ///
    /// compute residual to check solutions
    ///
    const double res = computeRelativeResidual(A, x, b);    

    std::cout << "For small matrices (test.mtx), the residual can be large; try with a bigger matrix\n";
    std::cout << "CuSparseTriSolve: residual = " << res << "\n\n";

  }
#else
  r_val = -1;
  std::cout << "CUDA is NOT configured in Trilinos" << std::endl;
#endif
  
  Kokkos::finalize();

  return r_val;
}
