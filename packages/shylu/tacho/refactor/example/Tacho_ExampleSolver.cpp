#include "Teuchos_CommandLineProcessor.hpp"

#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

int main (int argc, char *argv[]) {
  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Pardiso Chol algorithms on Kokkos::OpenMP execution space.\n");

  int nthreads = 1;
  clp.setOption("kokkos-threads", &nthreads, "Number of threads");

  int max_num_superblocks = 4;
  clp.setOption("max-num-superblocks", &max_num_superblocks, "Max number of superblocks");

  bool verbose = true;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");
  
  std::string file = "test.mtx";
  clp.setOption("file", &file, "Input file (MatrixMarket SPD matrix)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of RHS vectors");
  
  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  Kokkos::initialize(argc, argv);
  if (std::is_same<Kokkos::DefaultHostExecutionSpace,Kokkos::Serial>::value) 
    std::cout << "Kokkos::Serial\n";
  else
    Kokkos::DefaultHostExecutionSpace::print_configuration(std::cout, false);
  
  int r_val = 0;
  
  {
    /// basic typedef
    typedef int ordinal_type;
    //typedef size_t size_type; // not used here
    typedef double value_type;
    
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
    solver.setVerbose(verbose);
    solver.setMaxNumberOfSuperblocks(max_num_superblocks);

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

    std::cout << "TachoSolver: residual = " << res << std::endl;
  }
  Kokkos::finalize();

  return r_val;
}
