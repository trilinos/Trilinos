#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

#include "Tacho_CommandLineParser.hpp" 

using ordinal_type = Tacho::ordinal_type;

/// select a kokkos task scheudler
/// - DeprecatedTaskScheduler, DeprecatedTaskSchedulerMultiple
/// - TaskScheduler, TaskSchedulerMultiple, ChaseLevTaskScheduler
#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskScheduler<T>;
static const char * scheduler_name = "DeprecatedTaskScheduler";
#endif
#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskSchedulerMultiple<T>;
static const char * scheduler_name = "DeprecatedTaskSchedulerMultiple";
#endif
#if defined(TACHO_USE_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::TaskScheduler<T>;
static const char * scheduler_name = "TaskScheduler";
#endif
#if defined(TACHO_USE_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::TaskSchedulerMultiple<T>;
static const char * scheduler_name = "TaskSchedulerMultiple";
#endif
#if defined(TACHO_USE_CHASELEV_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::ChaseLevTaskScheduler<T>;
static const char * scheduler_name = "ChaseLevTaskScheduler";
#endif

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

  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

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

  const bool detail = false;

  typedef typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::device_type device_type;
  typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;

  typedef TaskSchedulerType<typename device_type::execution_space> scheduler_type;
  
  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace",   detail);
  printf("Scheduler Type = %s\n", scheduler_name);

  int r_val = 0;
  
  {
    /// crs matrix format and dense multi vector
    //typedef Tacho::CrsMatrixBase<value_type,device_type> CrsMatrixBaseType;
    typedef Tacho::CrsMatrixBase<value_type,host_device_type> CrsMatrixBaseTypeHost;
    
    typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> DenseMultiVectorType;
    //typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,host_device_type> DenseMultiVectorTypeHost;

    /// read a spd matrix of matrix market format
    CrsMatrixBaseTypeHost A;
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
    
    ///
    /// * to wrap triple pointers, declare following view types
    ///   typedef Kokkos::View<ordinal_type*,host_device_type>  ordinal_type_array;
    ///   typedef Kokkos::View<size_type*,   host_device_type>  size_type_array;
    ///   typedef Kokkos::View<value_type*,  host_device_type> value_type_array;
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
    Tacho::Solver<value_type,scheduler_type> solver;
    solver.setMatrixType(sym, posdef);
    solver.setVerbose(verbose);
    solver.setMaxNumberOfSuperblocks(max_num_superblocks);
    solver.setSmallProblemThresholdsize(small_problem_thres);
    solver.setBlocksize(mb);
    solver.setPanelsize(nb);

    auto values_on_device = Kokkos::create_mirror_view(typename device_type::memory_space(), A.Values());
    Kokkos::deep_copy(values_on_device, A.Values());

    /// inputs are used for graph reordering and analysis
    solver.analyze(A.NumRows(),
                   A.RowPtr(),
                   A.Cols());
    
    /// symbolic structure can be reused
    solver.factorize(values_on_device);
    
    DenseMultiVectorType 
      b("b", A.NumRows(), nrhs), // rhs multivector
      x("x", A.NumRows(), nrhs), // solution multivector
      t("t", A.NumRows(), nrhs); // temp workspace (store permuted rhs)
    
    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);
      Kokkos::fill_random(b, random, value_type(1));
    }
    solver.solve(x, b, t);
    
    const double res = solver.computeRelativeResidual(values_on_device, x, b);

    std::cout << "TachoSolver: residual = " << res << "\n\n";

    solver.release();
  }
  Kokkos::finalize();

  return r_val;
}
