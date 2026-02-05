// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#include "Kokkos_Random.hpp"

#include "Tacho_CommandLineParser.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"
#include "Tacho_MatrixMarket.hpp"

#ifdef TACHO_HAVE_TEUCHOS
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#endif

using ordinal_type = Tacho::ordinal_type;

template <typename value_type> int driver(int argc, char *argv[]) {
  using arith_traits = Tacho::ArithTraits<value_type>;
  using mag_type = typename arith_traits::mag_type;
#ifdef TACHO_HAVE_TEUCHOS
  using Teuchos::TimeMonitor;
  using Teuchos::StackedTimer;
  Teuchos::RCP<Teuchos::StackedTimer> stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Tacho_ExampleDriver"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
#endif

  int nthreads = 1;
  bool verbose = false;
  bool sanitize = false;
  bool duplicate = false;
  std::string file = "test.mtx";
  std::string rhs_file = "";
  std::string graph_file = "";
  std::string weight_file = "";
  int mm_base = 1;
  bool default_setup = false;
  int dofs_per_node = 1;
  int graph_algo_type = -1;
  bool order_connected_graph_separately = true;
#if defined(KOKKOS_ENABLE_HIP)
  bool storeTranspose = true;
#else
  bool storeTranspose = false;
#endif
  bool perturbPivot = false;
  int nrhs = 1;
  bool randomRHS = false;
  bool onesRHS = false;
  std::string method_name = "chol";
  int method = 1; // 0 - LDL no pivot, 1 - Chol, 2 - LDL, 3 - SymLU
  int small_problem_thres = 1024;
  int device_factor_thres = 64;
  int device_solve_thres = 128;
  int variant = 0;
  int nstreams = 8;
  bool no_warmup = false;
  int niters = 1;
  int nfacts = 2;
  int nsolves = 10;

  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);
  opts.set_option<bool>("duplicate", "Flag to duplicate input graph in the solver", &duplicate);
  opts.set_option<int>("mm-base", "Base for reading matrix file", &mm_base);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<std::string>("rhs", "Input RHS file", &rhs_file);
  opts.set_option<std::string>("graph", "Input condensed graph", &graph_file);
  opts.set_option<std::string>("weight", "Input condensed graph weight", &weight_file);
  opts.set_option<bool>("default", "Flag to use default parameters", &default_setup);
  opts.set_option<int>("dofs-per-node", "# DoFs per node", &dofs_per_node);
  opts.set_option<bool>("order-connected-graph", "Flag to order connected graph separately (METIS)", &order_connected_graph_separately);
  opts.set_option<int>("graph-algo", "Type of graph algorithm (0: Natural, 1: AMD, 2: METIS)", &graph_algo_type);
  opts.set_option<bool>("store-trans", "Flag to store transpose", &storeTranspose);
  opts.set_option<bool>("perturb", "Flag to perturb tiny pivots", &perturbPivot);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);
  opts.set_option<std::string>("method", "Solution method: ldl-nopiv, chol, ldl, lu", &method_name);
  opts.set_option<int>("small-problem-thres", "LAPACK is used smaller than this thres", &small_problem_thres);
  opts.set_option<int>("device-factor-thres", "Device function is used above this subproblem size",
                       &device_factor_thres);
  opts.set_option<int>("device-solve-thres", "Device function is used above this subproblem size", &device_solve_thres);
  opts.set_option<int>("variant", "algorithm variant in levelset scheduling; 0, 1 and 2", &variant);
  opts.set_option<int>("nstreams", "# of streams used in CUDA; on host, it is ignored", &nstreams);
  opts.set_option<bool>("one-rhs", "Set RHS to be ones", &onesRHS);
  opts.set_option<bool>("random-rhs", "Set RHS to be random", &randomRHS);
  opts.set_option<bool>("no-warmup", "Flag to turn off warmup", &no_warmup);
  opts.set_option<int>("niters", "# of iterations", &niters);
  opts.set_option<int>("nfacts", "# of factorizations to perform", &nfacts);
  opts.set_option<int>("nsolves", "# of solves to perform", &nsolves);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  if (method_name == "ldl-nopiv")
    method = 0;
  else if (method_name == "chol")
    method = 1;
  else if (method_name == "ldl")
    method = 2;
  else if (method_name == "lu")
    method = 3;
  else {
    std::cout << "Error: not supported solution method\n";
    return -1;
  }

  int r_val = 0;
  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;
    using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
    using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
    using CrsMatrixBaseTypeHost = Tacho::CrsMatrixBase<value_type, host_device_type>;
    using DenseMultiVectorType = Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type>;
    using size_type_array_host = typename CrsMatrixBaseTypeHost::size_type_array;
    using ordinal_type_array_host = typename CrsMatrixBaseTypeHost::ordinal_type_array;

    std::cout << std::endl << "   ---------------------- " << std::endl;
    Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
    Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace", detail);
    std::cout << "     Method Name:: " << method_name << std::endl;
    std::cout << "     Solver Type:: " << variant << std::endl;
    std::cout << "          # RHSs:: " << nrhs << std::endl << std::endl;
    if (default_setup) {
      std::cout << "   Using default Parameters " << std::endl;
    } else {
      std::cout << "   Using non default Parameters " << std::endl;
      std::cout << "       # Streams:: " << nstreams << std::endl;
      std::cout << "    Small Poblem:: " << small_problem_thres << std::endl;
      std::cout << "   Device Thresh:: " << device_factor_thres << ", " << device_solve_thres << std::endl;
    }
    std::cout << "  ---------------------- " << std::endl << std::endl;

    /// crs matrix format and dense multi vector
    /// read a spd matrix of matrix market format
    CrsMatrixBaseTypeHost A;
    {
      if (Tacho::MatrixMarket<value_type>::read(file, A, mm_base, sanitize, verbose) != 0) {
        return -1;
      }
      if (verbose) A.showMe(std::cout, false);
    }

    /// read graph file if available
    ordinal_type m_graph(0);
    size_type_array_host ap_graph;
    ordinal_type_array_host aw_graph, aj_graph;
    {
      if (graph_file.length() > 0 && weight_file.length() > 0) {
        {
          std::ifstream in;
          in.open(graph_file);
          if (in.good()) {
            std::cout << " > Condensed graph file: " << graph_file << std::endl;
          } else {
            std::cout << "Failed to open the graph file: " << graph_file << std::endl;
            return -1;
          }
          in >> m_graph;

          ap_graph = size_type_array_host("ap", m_graph + 1);
          for (ordinal_type i = 0, iend = m_graph + 1; i < iend; ++i)
            in >> ap_graph(i);

          aj_graph = ordinal_type_array_host("aj", ap_graph(m_graph));
          for (ordinal_type i = 0; i < m_graph; ++i) {
            const ordinal_type jbeg = ap_graph(i), jend = ap_graph(i + 1);
            for (ordinal_type j = jbeg; j < jend; ++j) {
              in >> aj_graph(j);
              aj_graph(j) --; // base-one
            }
          }
        }

        {
          std::ifstream in;
          in.open(weight_file);
          if (in.good()) {
            std::cout << " > Weight file for condensed graph: " << weight_file << std::endl;
          } else {
            std::cout << "Failed in open the file: " << weight_file << std::endl;
            return -1;
          }
          ordinal_type m(0);
          in >> m;
          in >> m_graph;
          aw_graph = ordinal_type_array_host("aw", m_graph);
          for (ordinal_type i = 0; i < m_graph; ++i)
            in >> aw_graph(i);
        }
      }
    }

    /// solver object
    Tacho::Driver<value_type, device_type> solver;

    /// common options
    solver.setVerbose(verbose);
    solver.setSolutionMethod(method);
    solver.setLevelSetOptionAlgorithmVariant(variant);

    if (!default_setup) {
      /// graph options
      if (order_connected_graph_separately) {
        solver.setOrderConnectedGraphSeparately();
      }
      if (graph_algo_type >= 0) {
        solver.setGraphAlgorithmType(graph_algo_type);
      }
      /// levelset options
      solver.setLevelSetOptionNumStreams(nstreams);
      solver.setSmallProblemThresholdsize(small_problem_thres);
      solver.setLevelSetOptionDeviceFunctionThreshold(device_factor_thres, device_solve_thres);
      if (perturbPivot) {
        if (verbose) std::cout << " > perturb tiny pivots" << std::endl;
        solver.useDefaultPivotTolerance();
      }
      solver.storeExplicitTranspose(storeTranspose);
      if (verbose) {
        if (storeTranspose) {
          std::cout << " > store transpose " << std::endl;
        }
        std::cout << std::endl;
      }
    }

    auto values_on_device = Kokkos::create_mirror_view(typename device_type::memory_space(), A.Values());
    Kokkos::deep_copy(values_on_device, A.Values());

    /// inputs are used for graph reordering and analysis
    {
#ifdef TACHO_HAVE_TEUCHOS
      Teuchos::TimeMonitor localTimer(*Teuchos::TimeMonitor::getNewTimer("Analyze"));
#endif
      if (m_graph > 0 && m_graph < A.NumRows())
        solver.analyze(A.NumRows(), A.RowPtr(), A.Cols(), m_graph, ap_graph, aj_graph, aw_graph);
      else if (dofs_per_node > 1) {
        if (verbose) std::cout << " > DoFs / node = " << dofs_per_node << std::endl;
        solver.analyze(A.NumRows(), dofs_per_node, A.RowPtr(), A.Cols());
      } else
        solver.analyze(A.NumRows(), A.RowPtr(), A.Cols());
    }

    /// create numeric tools and levelset tools
    Kokkos::Timer timer;
#ifdef TACHO_HAVE_TEUCHOS
    {
      Teuchos::TimeMonitor localTimer(*Teuchos::TimeMonitor::getNewTimer("Initialize"));
      solver.initialize();
    }
#else
    solver.initialize();
    double initi_time = timer.seconds();
#endif

    // setup vectors
    DenseMultiVectorType
        b("b", A.NumRows(), nrhs),  // rhs multivector
        x("x", A.NumRows(), nrhs),  // solution multivector
        t("t", A.NumRows(), nrhs);  // temp workspace (store permuted rhs)
    {
      if (rhs_file.length() > 0) {
        if(Tacho::MatrixMarket<value_type>::readDenseVectors(rhs_file, b, verbose) != 0) {
          return -1;
        }
      } else if (onesRHS) {
        if (verbose) std::cout << std::endl << " > RHS = ones" << std::endl << std::endl;
        const value_type one(1.0);
        Kokkos::deep_copy (b, one);
      } else if (randomRHS) {
        if (verbose) std::cout << std::endl << " > RHS = rands" << std::endl << std::endl;
        Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);
        Kokkos::fill_random(b, random, value_type(1));
      } else {
        if (verbose) std::cout << std::endl << " > RHS = A*ones" << std::endl << std::endl;
        const value_type one(1.0);
        Kokkos::deep_copy (x, one);
        solver.computeSpMV(values_on_device, x, b);
      }
    }

    bool success = true;
    mag_type tol = sqrt(arith_traits::epsilon()); // loose accuracy tol..
    for (int iter = 0; iter < niters; iter++) {
      try {
        /// symbolic structure can be reused
        std::cout << "\n === Iteration " << iter << " === " << "\n";

        /// factorize
        if (!no_warmup) {
          // warm-up
          solver.factorize(values_on_device);
        }

#ifdef TACHO_HAVE_TEUCHOS
        for (int i = 0; i < nfacts; ++i) {
          Teuchos::TimeMonitor localTimer(*Teuchos::TimeMonitor::getNewTimer("Factorize"));
          solver.factorize(values_on_device);
        }
#else
        timer.reset();
        for (int i = 0; i < nfacts; ++i) {
          solver.factorize(values_on_device);
        }
        double facto_time = timer.seconds();
#endif

        /// solve
        bool pass = true;
        double solve_time = 0.0;
        if (!no_warmup) {
          // warm-up
          timer.reset();
          solver.solve(x, b, t);
          solve_time = timer.seconds();
          const double res = solver.computeRelativeResidual(values_on_device, x, b);
          std::cout << "TachoSolver (warm-up): residual = " << res << " time " << solve_time << "\n";
          if (res > tol) pass = false;
        }

        for (int i = 0; i < nsolves; ++i) {
#ifdef TACHO_HAVE_TEUCHOS
          {
            Teuchos::TimeMonitor localTimer(*Teuchos::TimeMonitor::getNewTimer("Solve"));
            solver.solve(x, b, t);
          }
#else
          timer.reset();
          solver.solve(x, b, t);
          solve_time += timer.seconds();
#endif
          const mag_type res = solver.computeRelativeResidual(values_on_device, x, b);
          if (res > tol) pass = false;
          std::cout << "TachoSolver: residual = " << res << "\n";
        }
        if (!pass) success = false;
        std::cout << (pass ? " PASS" : " FAIL") << std::endl;
        std::cout << std::endl;
      } catch (const std::exception &e) {
        success = false;
        std::cerr << "Error: exception is caught: \n" << e.what() << "\n";
        std::cout << " EXCEPTION" << std::endl;
      }
    }
    if (verbose) {
      solver.printParameters();
    }
#ifdef TACHO_HAVE_TEUCHOS
    stackedTimer->stopBaseTimer();
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::rcp(new Teuchos::SerialComm<int>());
    // print stacked timer
    {
      Teuchos::StackedTimer::OutputOptions options;
      options.num_histogram=3;
      options.print_warnings = false;
      options.output_histogram = true;
      options.output_fraction=true;
      options.output_minmax = true;
      stackedTimer->report(std::cout, comm, options);
    }
    // generate watcher report
    if (!success) {
      std::cerr << "\n Error: Some of the residual norms were too large\n\n";
      stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Tacho_ExampleDriver"));
    }
    std::string testName = "Tacho_ExampleDriver_" + method_name + "_" + std::to_string(variant);
    auto xmlOut = stackedTimer->reportWatchrXML(testName, comm);
    if(xmlOut.length()) {
      std::cout << "\nAlso created Watchr performance report (" << testName << ") in " << xmlOut << '\n';
    } else {
      std::cout << "\nFailed to create Watchr performance report (" << testName << ") in " << xmlOut << '\n';
    }
    std::cout << std::endl;
#else
    std::cout << " Initi Time " << initi_time << std::endl;
    std::cout << " Facto Time " << facto_time / (double)nfacts << std::endl;
    std::cout << " Solve Time " << solve_time / (double)nsolves << std::endl;
#endif
    if (verbose) {
      std::cout << std::endl;
      std::cout << " > nnz = " << solver.getNumNonZerosU() << std::endl << std::endl;
    }
    std::cout << std::endl;
    solver.release();
  }
  Kokkos::finalize();

  return r_val;
}
