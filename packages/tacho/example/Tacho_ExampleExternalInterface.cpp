#include "Tacho_ExampleExternalInterface.hpp"
#include "Tacho_CommandLineParser.hpp" 

#if defined(TACHO_USE_INT_INT) 

double solutionError(const int numRows,
                     const double* rhs,
                     const double* sol,
                     const int* rowBegin,
                     const int* columns,
                     const double* values)
{
  double normRhsSquared(0), normErrorSquared(0);
  for (int i=0; i<numRows; i++) {
    double resid = rhs[i];
    for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
      const int col = columns[j];
      resid -= values[j]*sol[col];
    }
    double absVal = std::abs(rhs[i]);
    normRhsSquared += absVal*absVal;
    absVal = std::abs(resid);
    normErrorSquared += absVal*absVal;
  }

  double normError = std::sqrt(normErrorSquared);
  double normRhs = std::sqrt(normRhsSquared);
  return normError/normRhs;
}

void testTachoSolver(int numRows,
                     int* rowBegin,
                     int* columns,
                     double* values,
                     const int numRuns = 1)
{
  Kokkos::Impl::Timer timer;
  const double tol = 1e-10;//  tol *= 10; // original tolerance may be too tight

  ///
  /// Tacho options
  ///
  std::vector<int> tachoParams(tacho::INDEX_LENGTH);
  tachoParams[tacho::USEDEFAULTSOLVERPARAMETERS] = 0;
  tachoParams[tacho::VERBOSITY] = 0;
  tachoParams[tacho::SMALLPROBLEMTHRESHOLDSIZE] = 1024;
  
#if defined (KOKKOS_ENABLE_CUDA)
  tachoParams[tacho::TASKING_OPTION_MAXNUMSUPERBLOCKS] = 32;
  tachoParams[tacho::TASKING_OPTION_BLOCKSIZE] = 64;
  tachoParams[tacho::TASKING_OPTION_PANELSIZE] = 32;

  tachoParams[tacho::LEVELSET_OPTION_SCHEDULING] = 1;
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_LEVEL_CUT] = 0;
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_FACTOR_THRES] = 64;
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_SOLVE_THRES] = 128;
  tachoParams[tacho::LEVELSET_OPTION_NSTREAMS] = 8;
#else
#  ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  tachoParams[tacho::TASKING_OPTION_MAXNUMSUPERBLOCKS] = std::max(tacho::host_space::thread_pool_size(0)/2,1);
#  else
  tachoParams[tacho::TASKING_OPTION_MAXNUMSUPERBLOCKS] = std::max(tacho::host_space::impl_thread_pool_size(0)/2,1);
#  endif
  tachoParams[tacho::TASKING_OPTION_BLOCKSIZE] = 256;
  tachoParams[tacho::TASKING_OPTION_PANELSIZE] = 128;

  tachoParams[tacho::LEVELSET_OPTION_SCHEDULING] = 0;
  // the following options are not used and set dummy values
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_LEVEL_CUT] = 0;
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_FACTOR_THRES] = 0;
  tachoParams[tacho::LEVELSET_OPTION_DEVICE_SOLVE_THRES] = 0;
  tachoParams[tacho::LEVELSET_OPTION_NSTREAMS] = 0;
#endif

  ///
  /// Tacho solver analyze + factorize (fence is required)
  ///
  tacho::tachoSolver<double> solver(tachoParams.data());
  solver.Initialize(numRows, rowBegin, columns, values);
  Kokkos::fence();


  ///
  /// Export supernodes
  ///
  std::vector<int> supernodes;
  solver.exportSupernodes(supernodes);

  ///
  /// Export matrix and permutation vector
  ///
  std::vector<int> rowBeginU;
  std::vector<int> columnsU;
  std::vector<int> perm;
  std::vector<double> valuesU;

  solver.exportUpperTriangularFactorsToCrsMatrix(rowBeginU,
                                                 columnsU,
                                                 valuesU,
                                                 perm);
  ///
  /// std vector right hand side
  /// if an application uses std vector for interfacing rhs,
  /// it requires additional copy. it is better to directly 
  /// use a kokkos device view.
  ///
  std::vector<double> rhs(numRows), sol(numRows);
  { /// randomize rhs
    const unsigned int seed = 0;
    srand(seed);
    for (int i=0;i<numRows;++i)
      rhs[i] = rand()/((double)RAND_MAX+1.0);
  }
  
  /// this example only works for single right hand side
  const int NRHS = 1;
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, tacho::device_type> ViewVectorType;
  ViewVectorType x("x", numRows, NRHS);

#if defined (KOKKOS_ENABLE_CUDA)
  /// transfer b into device
  ViewVectorType b(Kokkos::ViewAllocateWithoutInitializing("b"), numRows, NRHS);
  Kokkos::deep_copy(Kokkos::subview(b, Kokkos::ALL(), 0), 
                    Kokkos::View<double*,tacho::device_type>(rhs.data(), numRows));
#else
  /// wrap rhs data with view
  ViewVectorType b(rhs.data(), numRows, NRHS);
#endif

  timer.reset();
  
  for (int run=0; run<numRuns; run++) {
    solver.MySolve(NRHS, b, x);
    Kokkos::fence();
  }
  double runTime = timer.seconds();

  /// computing residual on host
  auto h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  
  std::cout << "ExternalInterface:: solve time (" << numRuns << ") = " << runTime << std::endl;
  double relativeError = solutionError(numRows, rhs.data(), h_x.data(), rowBegin, columns, values);
  std::cout << "ExternalInterface:: relative error = " << relativeError << std::endl;
  if (relativeError > tol) {
    std::cout << "ExternalInterface:: tacho FAILS (rel Error = " 
              << relativeError << ", tol = " << tol << std::endl;
  } else {
    std::cout << "ExternalInterface:: tacho PASS " << std::endl;
  }
}

int main(int argc, char *argv[]) {
  Tacho::CommandLineParser opts("This example shows potential use case of Tacho interface");

  std::string file = "test.mtx";
  int niter = 1;
  bool verbose = true;
  bool sanitize = false;

  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("niter", "# of solver iterations", &niter);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);
  const bool detail = false;
  Tacho::printExecSpaceConfiguration<tacho::exec_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<tacho::host_space>("HostSpace",   detail);
  {

    Tacho::CrsMatrixBase<double,tacho::host_device_type> A;
    {
      std::ifstream in;
      in.open(file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file << std::endl;
        Kokkos::finalize();
        return -1;
      }
    }
    Tacho::MatrixMarket<double>::read(file, A, sanitize, verbose);

    int numRows = A.NumRows(), *rowBegin = A.RowPtr().data(), *columns = A.Cols().data();
    double *values = A.Values().data();

    testTachoSolver(numRows, rowBegin, columns, values, niter);
  }
  Kokkos::finalize();
}

#else
#include <iostream>
int main(int argc, char *argv[]) {
  std::cout << "ExternalInterface example is NOT enabled; please configure with -D Tacho_ENABLE_INT_INT=ON\n";
  return 0;
}
#endif
