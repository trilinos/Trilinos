// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_set>

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include <KokkosKernels_TestStringUtils.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include "KokkosSparse_IOUtils.hpp"
#include <KokkosKernels_config.h>
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_TestMatrixUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_SortCrs.hpp"

using std::cout;
using std::string;
using namespace KokkosSparse;

static char* getNextArg(int& i, int argc, char** argv) {
  i++;
  if (i >= argc) {
    std::cerr << "Error: expected additional command-line argument!\n";
    exit(1);
  }
  return argv[i];
}

struct Sort_Parameters {
  const char* matrix_path = nullptr;
  int repeat              = 10;
  int n                   = 10000;
  int nnzPerRow           = 27;
  int numLongRows         = 0;
  int minNnzPerLongRow    = 1000;
  int maxNnzPerLongRow    = 2000;
  bool alreadySorted      = false;
  // Whether to shuffle entries of a matrix read from file.
  // Entries of random generated matrix are always shuffled.
  bool shuffle       = false;
  bool graphOnly     = false;
  SortAlgorithm algo = SortAlgorithm::DEFAULT;
};

template <typename crsMat_t>
crsMat_t generateLongRowMatrix(const Sort_Parameters& params, typename crsMat_t::ordinal_type ncols) {
  typedef typename crsMat_t::value_type scalar_t;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::size_type size_type;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::index_type::non_const_type entries_view_t;
  typedef typename crsMat_t::row_map_type::non_const_type rowmap_view_t;
  // Generate random diag. dominant matrix
  srand(245);
  std::vector<size_type> rowmap = {0};
  std::vector<lno_t> entries;
  std::vector<scalar_t> values;
  std::vector<lno_t> rowLengths;
  lno_t numRows = params.n;
  for (lno_t i = 0; i < numRows; i++) {
    if (i < params.numLongRows) {
      lno_t interval = params.maxNnzPerLongRow - params.minNnzPerLongRow;
      lno_t rowLen;
      if (interval == 0)
        rowLen = params.maxNnzPerLongRow;
      else
        rowLen = params.minNnzPerLongRow + rand() % interval;
      if (rowLen > numRows) rowLen = numRows;
      rowLengths.push_back(rowLen);
    } else
      rowLengths.push_back(params.nnzPerRow);
  }
  std::shuffle(rowLengths.begin(), rowLengths.end(), std::mt19937(std::random_device()()));
  size_type totalEntries = 0;
  int randSteps          = 1000000;
  // Set of columns inserted so far into current short row
  std::unordered_set<lno_t> shortRowEntries;
  // Set of all possible rows, randomly permuted (select long row entries from
  // head)
  std::vector<lno_t> longRowEntries(numRows);
  for (lno_t i = 0; i < numRows; i++) longRowEntries[i] = i;
  const scalar_t one = 1;
  for (lno_t i = 0; i < numRows; i++) {
    shortRowEntries.clear();
    bool rowIsLong = rowLengths[i] > params.nnzPerRow;
    if (rowIsLong) std::shuffle(longRowEntries.begin(), longRowEntries.end(), std::mt19937(std::random_device()()));
    for (lno_t ent = 0; ent < rowLengths[i]; ent++) {
      if (ent == 0) {
        entries.push_back(i);
        values.push_back(5.0 + 3.0 * (rand() % randSteps) / randSteps * one);
      } else {
        if (rowIsLong)
          entries.push_back(longRowEntries[ent]);
        else {
          // re-roll random column until one is found that isn't already in row
          lno_t col;
          while (true) {
            col = rand() % ncols;
            if (shortRowEntries.find(col) == shortRowEntries.end()) {
              shortRowEntries.insert(col);
              break;
            }
          }
          entries.push_back(col);
        }
        values.push_back((-0.1 + (0.2 * (rand() % randSteps) / randSteps)) * one);
      }
    }
    totalEntries += rowLengths[i];
    rowmap.push_back(totalEntries);
  }
  scalar_view_t valuesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values"), totalEntries);
  entries_view_t entriesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries"), totalEntries);
  rowmap_view_t rowmapView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Rowmap"), numRows + 1);
  Kokkos::deep_copy(valuesView, Kokkos::View<scalar_t*, Kokkos::HostSpace>(values.data(), totalEntries));
  Kokkos::deep_copy(entriesView, Kokkos::View<lno_t*, Kokkos::HostSpace>(entries.data(), totalEntries));
  Kokkos::deep_copy(rowmapView, Kokkos::View<size_type*, Kokkos::HostSpace>(rowmap.data(), numRows + 1));
  crsMat_t A("A", numRows, numRows, totalEntries, valuesView, rowmapView, entriesView);
  A = KokkosSparse::sort_and_merge_matrix(A);
  return A;
}

template <typename device_t>
void runSort(const Sort_Parameters& params) {
  using scalar_t  = KokkosKernels::default_scalar;
  using lno_t     = KokkosKernels::default_lno_t;
  using size_type = KokkosKernels::default_size_type;
  lno_t ncols     = params.n;
  if (!params.matrix_path) {
    // For generating matrix, params.n is number of rows.
    // By default the matrix is square, but if needed adjust number of columns
    // to fit requested nz per row if n is too small.
    lno_t minCols = 2 * (params.numLongRows ? params.maxNnzPerLongRow : params.nnzPerRow);

    if (ncols < minCols) {
      std::cout << "Increasing number of columns to " << minCols << " to fit requested nz/row.\n";
      ncols = minCols;
    }
  }
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type> crsMat_t;
  crsMat_t A;
  bool doShuffle = params.shuffle || (params.matrix_path == nullptr);
  bool doPreSort = params.alreadySorted;
  if (params.matrix_path) {
    std::cout << "Reading matrix from file...";
    std::cout.flush();
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.matrix_path);
    std::cout << "done" << std::endl;
  } else {
    std::cout << "Generating matrix...";
    std::cout.flush();
    A = generateLongRowMatrix<crsMat_t>(params, ncols);
    std::cout << "done" << std::endl;
  }
  if (A.numRows() == 0 || A.numCols() == 0) {
    std::cout << "Matrix has 0 rows or 0 cols; sorting is a no-op\n";
    return;
  }
  // Copy of entries to be used for input each repetition (either sorted or shuffled in each row)
  typename crsMat_t::index_type inputEntries("input entries", A.nnz());
  // Randomly shuffle or sort indices within each row of Acopy, depending on params.alreadySorted
  // We don't care about what the output values are for testing performance, so don't bother to correspondingly shuffle
  // them.
  {
    int m        = A.numRows();
    auto rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
    auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
    std::random_device rd;
    std::mt19937 g(rd());
    for (int i = 0; i < m; i++) {
      size_t rowBegin = rowmap(i);
      size_t rowEnd   = rowmap(i + 1);
      if (doPreSort)
        std::sort(entries.data() + rowBegin, entries.data() + rowEnd);
      else if (doShuffle)
        std::shuffle(entries.data() + rowBegin, entries.data() + rowEnd, g);
    }
    Kokkos::deep_copy(inputEntries, entries);
  }
  Kokkos::Timer timer;
  double sortTime = 0;
  for (int i = 0; i < params.repeat; i++) {
    // Set up input matrix again.
    Kokkos::deep_copy(A.graph.entries, inputEntries);
    timer.reset();
    if (params.graphOnly) {
      KokkosSparse::sort_crs_graph(A.graph, A.numCols(), params.algo);
    } else {
      KokkosSparse::sort_crs_matrix(A, params.algo);
    }
    Kokkos::fence();
    sortTime += timer.seconds();
  }

  if (params.graphOnly)
    std::cout << "\n*** Mean graph sort time over " << params.repeat << " repeats: ";
  else
    std::cout << "\n*** Mean matrix sort time over " << params.repeat << " repeats: ";
  std::cout << sortTime / params.repeat << " s\n";
}

int main(int argc, char** argv) {
  // Expect two args: matrix name and device flag.
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    cout << "* General options:\n";
    cout << "--serial, --openmp, --cuda, --hip, --sycl, --threads: select backend\n";
    cout << "--repeat N: sort matrix/graph N times and average time (default=10)\n";
    cout << "--mtx <file>: provide MatrixMarket file to read matrix from\n\n";
    cout << "* Randomized matrix settings, if not reading from file:\n";
    cout << "  --n <N> : number of rows/columns\n";
    cout << "  --nnz <N> : number of nonzeros in each regular row\n";
    cout << "  --long-rows <N> : number of long rows\n";
    cout << "  --min-long-row-nnz <N> : min number of nonzeros in each long "
            "row (default 1000)\n";
    cout << "  --max-long-row-nnz <N> : max number of nonzeros in each long "
            "row (default 2000)\n\n";

    cout << "* Sorting options:\n";
    cout << "--graph-only: sort CRS graph only, not matrix (default off)\n";
    cout << "--shuffle: (when reading from file only) randomly shuffle entries (default off)\n";
    cout << "--already-sorted: test with input already being sorted (default off)\n";
    cout << "Sort algorithm:\n";
    cout << "  --algo <algo> (default: DEFAULT)\n";
    cout << "         DEFAULT: choose based on heuristics\n";
    cout << "         BULK_SORT: sort all rows at once; better for very long rows\n";
    cout << "         RADIX: serial radix in each row (CPU only)\n";
    cout << "         SHELL: serial Shell in each row (CPU only)\n";
    return 0;
  }
  Kokkos::initialize(argc, argv);
  // device is just the name of the execution space, lowercase
  string deviceName;
  Sort_Parameters params;
  int i = 1;
  for (; i < argc; i++) {
    if (!strcmp(argv[i], "--mtx"))
      params.matrix_path = getNextArg(i, argc, argv);
    else if (!strcmp(argv[i], "--n"))
      params.n = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--nnz"))
      params.nnzPerRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--long-rows"))
      params.numLongRows = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--min-long-row-nnz"))
      params.minNnzPerLongRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--max-long-row-nnz"))
      params.maxNnzPerLongRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--serial"))
      deviceName = "serial";
    else if (!strcmp(argv[i], "--openmp"))
      deviceName = "openmp";
    else if (!strcmp(argv[i], "--threads"))
      deviceName = "threads";
    else if (!strcmp(argv[i], "--cuda"))
      deviceName = "cuda";
    else if (!strcmp(argv[i], "--hip"))
      deviceName = "hip";
    else if (!strcmp(argv[i], "--sycl"))
      deviceName = "sycl";
    else if (!strcmp(argv[i], "--repeat"))
      params.repeat = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--graph-only"))
      params.graphOnly = true;
    else if (!strcmp(argv[i], "--already-sorted"))
      params.alreadySorted = true;
    else if (!strcmp(argv[i], "--shuffle"))
      params.shuffle = true;
    else if (!strcmp(argv[i], "--algo")) {
      const char* algo = getNextArg(i, argc, argv);
      if (!strcmp(algo, "DEFAULT"))
        params.algo = SortAlgorithm::DEFAULT;
      else if (!strcmp(algo, "BULK_SORT"))
        params.algo = SortAlgorithm::BULK_SORT;
      else if (!strcmp(algo, "RADIX"))
        params.algo = SortAlgorithm::RADIX;
      else if (!strcmp(algo, "SHELL"))
        params.algo = SortAlgorithm::SHELL;
      else {
        std::cout << "Error: invalid sort algorithm.\n";
        Kokkos::finalize();
        exit(1);
      }
    } else {
      cout << "Error: unknown argument " << argv[i] << '\n';
      Kokkos::finalize();
      exit(1);
    }
  }
  bool run = false;
  if (!deviceName.length()) {
    std::cout << "Will run on default backend (" << Kokkos::DefaultExecutionSpace().name() << ")\n";
    runSort<Kokkos::DefaultExecutionSpace>(params);
    run = true;
  } else {
    std::cout << "Will run on requested backend (" << deviceName << ")\n";
  }
#ifdef KOKKOS_ENABLE_SERIAL
  if (deviceName == "serial") {
    runSort<Kokkos::Serial>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  if (deviceName == "openmp") {
    runSort<Kokkos::OpenMP>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_THREADS
  if (deviceName == "threads") {
    runSort<Kokkos::Threads>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_CUDA
  if (deviceName == "cuda") {
    runSort<Kokkos::Cuda>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_HIP
  if (deviceName == "hip") {
    runSort<Kokkos::HIP>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_SYCL
  if (deviceName == "sycl") {
    runSort<Kokkos::SYCL>(params);
    run = true;
  }
#endif
  if (!run) {
    std::cerr << "Error: device " << deviceName << " was requested but it's not enabled in this build.\n";
    return 1;
  }
  Kokkos::finalize();
  return 0;
}
