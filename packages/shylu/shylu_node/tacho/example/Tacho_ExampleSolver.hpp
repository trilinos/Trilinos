#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

#include "TachoExp_CommandLineParser.hpp" 

using ordinal_type = Tacho::Experimental::ordinal_type;

template<typename value_type>
int driver (int argc, char *argv[]) {
  int nthreads = 1;
  int max_num_superblocks = 4;
  bool verbose = true;
  std::string file = "test.mtx";
  int nrhs = 1;
  int sym = 3;
  int posdef = 1;
  int small_problem_thres = 1024;
  int mb = -1;
  int nb = -1;

  Tacho::Experimental::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<int>("max-num-superblocks", "Max number of superblocks", &max_num_superblocks);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("symmetric", "Symmetric type: 0 - unsym, 1 - structure sym, 2 - symmetric, 3 - hermitian", &sym);
  opts.set_option<int>("posdef", "Positive definite: 0 - indef, 1 - positive definite", &posdef);
  opts.set_option<int>("small-problem-thres", "LAPACK is used smaller than this thres", &small_problem_thres);
  opts.set_option<int>("mb", "Internal block size", &mb);
  opts.set_option<int>("nb", "Internal panel size", &nb);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);
  if (std::is_same<Kokkos::DefaultHostExecutionSpace,Kokkos::Serial>::value) 
    std::cout << "Kokkos::Serial\n";
  else
    Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);
  
  int r_val = 0;
  
  {
    /// basic typedef
    //typedef int ordinal_type;
    //typedef size_t size_type; // not used here
    //typedef double value_type;
    
    /// crs matrix format and dense multi vector
    typedef Tacho::Experimental::CrsMatrixBase<value_type,Kokkos::DefaultHostExecutionSpace> CrsMatrixBaseType;
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,Kokkos::DefaultHostExecutionSpace> DenseMultiVectorType;

    /// read a spd matrix of matrix market format
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
    
    ///
    /// * to wrap triple pointers, declare following view types
    ///   typedef Kokkos::View<ordinal_type*,Kokkos::DefaultHostExecutionSpace>  ordinal_type_array;
    ///   typedef Kokkos::View<size_type*,   Kokkos::DefaultHostExecutionSpace>  size_type_array;
    ///   typedef Kokkos::View<value_type*   ,Kokkos::DefaultHostExecutionSpace> value_type_array;
    ///
    /// * or, these can be derived from CrsMatrixBaseType
    ///   typedef typename CrsMatrixBaseType::ordinal_type_array ordinal_type_array;
    ///   typedef typename CrsMatrixBaseType::size_type_array    size_type_array;
    ///   typedef typename CrsMatrixBaseType::value_type_array   value_type_array;
    ///
    /// * wrap triple pointers (row_ptr, colidx_ptr, value_ptr) with views 
    ///   size_type_array ap(row_ptr, nrows + 1);
    ///   ordinal_type_array aj(colidx_ptr, nnz);
    ///   value_type_array ax(value_ptr, nnz);
    /// 
    /// * attach views into csr matrix 
    ///   CrsMatrixBaseType A;
    ///   A.setExternalMatrix(nrows, ncols, nnzm ap, aj, ax);
    ///  
    Tacho::Solver<value_type,Kokkos::DefaultHostExecutionSpace> solver;
    solver.setMatrixType(sym, posdef);
    solver.setVerbose(verbose);
    solver.setMaxNumberOfSuperblocks(max_num_superblocks);
    solver.setSmallProblemThresholdsize(small_problem_thres);
    solver.setBlocksize(mb);
    solver.setPanelsize(nb);

    /// inputs are used for graph reordering and analysis
    solver.analyze(A.NumRows(),
                   A.RowPtr(),
                   A.Cols());
    
    /// symbolic structure can be reused
    solver.factorize(A.Values());

    DenseMultiVectorType 
      b("b", A.NumRows(), nrhs), // rhs multivector
      x("x", A.NumRows(), nrhs), // solution multivector
      t("t", A.NumRows(), nrhs); // temp workspace (store permuted rhs)

    {
      Tacho::Experimental::Random<value_type> random;
      const ordinal_type m = A.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs)
        for (ordinal_type i=0;i<m;++i) 
          b(i, rhs) = random.value();
    }
    
    solver.solve(x, b, t);
    
    const double res = solver.computeRelativeResidual(A.Values(), x, b);

    std::cout << "TachoSolver: residual = " << res << "\n\n";

    solver.release();
  }
  Kokkos::finalize();

  return r_val;
}
