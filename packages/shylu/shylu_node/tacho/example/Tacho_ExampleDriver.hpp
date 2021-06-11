#include "Kokkos_Random.hpp"

#include "Tacho_Driver.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"
#include "Tacho_CommandLineParser.hpp" 

using ordinal_type = Tacho::ordinal_type;

template<typename value_type>
int driver (int argc, char *argv[]) {
  int nthreads = 1;
  bool verbose = true;
  bool sanitize = false;
  bool duplicate = false;
  std::string file = "test.mtx";
  std::string graph_file = "";
  std::string weight_file = "";
  int nrhs = 1;
  int sym = 2;
  int posdef = 1;
  int small_problem_thres = 1024;
  int device_factor_thres = 64;
  int device_solve_thres = 128;
  int variant = 0;
  int nstreams = 8;

  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);
  opts.set_option<bool>("duplicate", "Flag to duplicate input graph in the solver", &duplicate);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<std::string>("graph", "Input condensed graph", &graph_file);
  opts.set_option<std::string>("weight", "Input condensed graph weight", &weight_file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<int>("symmetric", "Symmetric type: 0 - unsym, 1 - structure sym, 2 - symmetric", &sym);
  opts.set_option<int>("posdef", "Positive definite: 0 - indef, 1 - positive definite", &posdef);
  opts.set_option<int>("small-problem-thres", "LAPACK is used smaller than this thres", &small_problem_thres);
  opts.set_option<int>("device-factor-thres", "Device function is used above this subproblem size", &device_factor_thres);
  opts.set_option<int>("device-solve-thres", "Device function is used above this subproblem size", &device_solve_thres);
  opts.set_option<int>("variant", "algorithm variant in levelset scheduling; 0, 1 and 2", &variant);
  opts.set_option<int>("nstreams", "# of streams used in CUDA; on host, it is ignored", &nstreams);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);

  const bool detail = false;

  using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
  using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;

  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace",   detail);

  int r_val = 0;  
  {
    /// crs matrix format and dense multi vector
    using CrsMatrixBaseTypeHost = Tacho::CrsMatrixBase<value_type,host_device_type>;
    using DenseMultiVectorType = Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type>;

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
      Tacho::MatrixMarket<value_type>::read(file, A, sanitize, verbose);
    }

    /// read graph file if available
    using size_type_array_host = typename CrsMatrixBaseTypeHost::size_type_array;
    using ordinal_type_array_host = typename CrsMatrixBaseTypeHost::ordinal_type_array;
    
    ordinal_type m_graph(0);
    size_type_array_host ap_graph;
    ordinal_type_array_host aw_graph, aj_graph;
    {
      if (graph_file.length() > 0 && weight_file.length() > 0) {
        {
          std::ifstream in;
          in.open(graph_file);
          if (!in.good()) {
            std::cout << "Failed in open the file: " << graph_file << std::endl;
            return -1;
          }
          in >> m_graph;
          
          ap_graph = size_type_array_host("ap", m_graph+1);
          for (ordinal_type i=0,iend=m_graph+1;i<iend;++i)
            in >> ap_graph(i);
          
          aj_graph = ordinal_type_array_host("aj", ap_graph(m_graph));
          for (ordinal_type i=0;i<m_graph;++i) {
            const ordinal_type jbeg = ap_graph(i), jend = ap_graph(i+1);
            for (ordinal_type j=jbeg;j<jend;++j)
              in >> aj_graph(j);
          }
        }
        
        {
          std::ifstream in;
          in.open(weight_file);
          if (!in.good()) {
            std::cout << "Failed in open the file: " << weight_file << std::endl;
            return -1;
          }
          ordinal_type m(0);
          in >> m;
          in >> m_graph;
          aw_graph = ordinal_type_array_host("aw", m_graph);
          for (ordinal_type i=0;i<m_graph;++i) 
            in >> aw_graph(i);
        }
      }
    }
    
    Tacho::Driver<value_type,device_type> solver;

    /// common options
    solver.setMatrixType(sym, posdef);
    solver.setSmallProblemThresholdsize(small_problem_thres);
    solver.setVerbose(verbose);

    /// graph options
    solver.setOrderConnectedGraphSeparately();

    /// levelset options
    solver.setLevelSetOptionDeviceFunctionThreshold(device_factor_thres, device_solve_thres);
    solver.setLevelSetOptionAlgorithmVariant(variant);
    solver.setLevelSetOptionNumStreams(nstreams);

    auto values_on_device = Kokkos::create_mirror_view(typename device_type::memory_space(), A.Values());
    Kokkos::deep_copy(values_on_device, A.Values());

    /// inputs are used for graph reordering and analysis
    if (m_graph > 0 && m_graph < A.NumRows())
      solver.analyze(A.NumRows(),
                     A.RowPtr(),
                     A.Cols(),
                     m_graph,
                     ap_graph,
                     aj_graph,
                     aw_graph);
    else 
      solver.analyze(A.NumRows(),
                     A.RowPtr(),
                     A.Cols());

    /// create numeric tools and levelset tools
    solver.initialize();

    /// symbolic structure can be reused
    for (int i=0;i<2;++i)
      solver.factorize(values_on_device);
    
    DenseMultiVectorType 
      b("b", A.NumRows(), nrhs), // rhs multivector
      x("x", A.NumRows(), nrhs), // solution multivector
      t("t", A.NumRows(), nrhs); // temp workspace (store permuted rhs)
    
    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);
      Kokkos::fill_random(b, random, value_type(1));
    }

    for (int i=0;i<3;++i)
      solver.solve(x, b, t);
    
    const double res = solver.computeRelativeResidual(values_on_device, x, b);

    std::cout << "TachoSolver: residual = " << res << "\n\n";

    solver.release();
  }
  Kokkos::finalize();

  return r_val;
}
