// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <cmath>
#include <numeric>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

using SC = Tpetra::CrsMatrix<>::scalar_type;
using GST = Tpetra::global_size_t;
using map_type = Tpetra::Map<>;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;
using crs_graph_type = Tpetra::CrsGraph<>;
using crs_matrix_type = Tpetra::CrsMatrix<SC>;
using export_type = Tpetra::Export<>;

struct CmdLineArgs {
  double densityFraction = 1.0;
  int lclNumRows = 10;
  int lclNumCols = -1;
  int numTrials = 100;
  int numOverlapRows = 3;
  bool contiguousMaps = true;
};

std::unique_ptr<GO[]>
getTargetRowMapIndices(const LO lclNumRows,
                       const int myRank,
                       const int numProcs)
{
  const GO gblNumRows = GO(lclNumRows) * GO(numProcs);
  const GO indexBase = 0;
  const GO gidBase = indexBase + GO(lclNumRows) * GO(myRank);

  // LIDs [0, N-1] -> GIDs gidBase + [N-1, 1, N-2, 2, ...].
  std::unique_ptr<GO[]> tgtGids(new GO[lclNumRows]);
  for (LO lid = 0; lid < lclNumRows; ++lid) {
    const LO offset = (lid % LO(2) == 0) ?
      lclNumRows - LO(1) - lid : lid;
    const GO gblRow = gidBase + GO(offset);
    TEUCHOS_ASSERT(gblRow >= indexBase);
    TEUCHOS_ASSERT(gblRow < indexBase + gblNumRows);
    tgtGids[lid] = gblRow;
  }
  // The original return using std::move (commented out below) returns the
  // following warning with gcc 9.2.0:
  // waring: moving a local object in a return statement prevents copy elision [-Wpessimizing-move]
  //return std::move(tgtGids);
  return tgtGids;
}

RCP<const map_type>
getTargetRowMap(const RCP<const Teuchos::Comm<int>>& comm,
                const CmdLineArgs& args)
{
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  const LO lclNumRows = args.lclNumRows;
  const GO gblNumRows = GO(lclNumRows) * GO(numProcs);
  const GO indexBase = 0;

  if(args.contiguousMaps) {
    return rcp(new map_type(gblNumRows, indexBase, comm));
  }
  else {
    std::unique_ptr<GO[]> tgtGids =
      getTargetRowMapIndices(lclNumRows, myRank, numProcs);
    return rcp(new map_type(gblNumRows, tgtGids.get(),
                            lclNumRows, indexBase, comm));
  }
}

// Only some of the target matrix's local rows get modified.  The
// Export object tells us which rows of the target matrix receive data
// from the source matrix.  I call these "incoming rows" in the code
// below.  The returned View is a View of bool, though I'm using int
// to avoid potential CUDA bool drama.  The View tells me, for a given
// local row of the target matrix, whether the doExport modified that
// local row.
Kokkos::View<const int*, crs_matrix_type::device_type>
getLclRowsToTest(const export_type& exporter)
{
  Teuchos::ArrayView<const LO> incomingRows = exporter.getRemoteLIDs();
  const map_type& tgtMap = *(exporter.getTargetMap());
  const LO tgtLclNumRows = LO(tgtMap.getLocalNumElements());

  using device_type = crs_matrix_type::device_type;
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  Kokkos::View<int*, device_type> lclRowsToTest
    (view_alloc("lclRowsToTest", WithoutInitializing), tgtLclNumRows);
  auto lclRowsToTest_h =
    Kokkos::create_mirror_view(Kokkos::HostSpace(), lclRowsToTest);
  Kokkos::deep_copy(lclRowsToTest_h, 0);

  for(const LO incomingRow : incomingRows) {
    lclRowsToTest_h[incomingRow] = 1;
  }

  const bool verbose = Tpetra::Details::Behavior::verbose();
  if(verbose) {
    using std::cerr;
    using std::endl;

    const auto comm = exporter.getSourceMap()->getComm();
    const int myRank = comm->getRank();
    if(myRank == 0) {
      cerr << "lclRowsToTest:" << endl;
    }
    std::ostringstream os;
    os << "Proc " << myRank << ": [";
    const LO numLclRowsToTest(lclRowsToTest_h.extent(0));
    for(LO k = 0; k < numLclRowsToTest; ++k) {
      os << lclRowsToTest_h[k];
      if (k + LO(1) < numLclRowsToTest) {
        os << ", ";
      }
    }
    os << "]" << endl;
    Tpetra::Details::gathervPrint(cerr, os.str(), *comm);
  }

  Kokkos::deep_copy(lclRowsToTest, lclRowsToTest_h);
  return lclRowsToTest;
}

RCP<const map_type>
getSourceRowMap(const RCP<const map_type>& tgtMap,
                const CmdLineArgs& args)
{
  const auto comm = tgtMap->getComm();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  const LO tgtLclNumRows = LO(args.lclNumRows);
  const GO tgtGblNumRows = GO(tgtLclNumRows) * GO(numProcs);
  const GO indexBase = tgtMap->getIndexBase();
  const GO tgtGidBase = indexBase + GO(tgtLclNumRows) * GO(myRank);
  const GO srcGidBase = indexBase +
    ((tgtGidBase + GO(tgtLclNumRows)) % tgtGblNumRows);
  const LO srcLclNumRows = LO(args.numOverlapRows);

  TEUCHOS_ASSERT(srcGidBase >= indexBase);
  TEUCHOS_ASSERT(srcGidBase < indexBase + tgtGblNumRows);

  // Construct Source Map so that copyAndPermute has nothing to do.
  // This should help focus the benchmark on pack and unpack.
  std::unique_ptr<GO[]> srcGids(new GO[srcLclNumRows]);
  for(LO lid = 0; lid < srcLclNumRows; ++lid) {
    const GO gblRow = srcGidBase + GO(lid);
    TEUCHOS_ASSERT(gblRow >= indexBase);
    TEUCHOS_ASSERT(gblRow < indexBase + tgtGblNumRows);
    srcGids[lid] = gblRow;
  }
  using Teuchos::OrdinalTraits;
  return rcp(new map_type(OrdinalTraits<GST>::invalid(), srcGids.get(),
                          srcLclNumRows, tgtMap->getIndexBase(), comm));
}

RCP<const map_type>
getDomainMap(const RCP<const Teuchos::Comm<int>>& comm,
             const CmdLineArgs& args)
{
  const int numProcs = comm->getSize();
  const GO gblNumCols = GO(args.lclNumCols) * GO(numProcs);
  const GO indexBase = 0;

  return rcp(new map_type(gblNumCols, indexBase, comm));
}

// Create a new timer with the given name if it hasn't already been
// created, else get the previously created timer with that name.
RCP<Time> getTimer(const std::string& timerName) {
  RCP<Time> timer = TimeMonitor::lookupCounter(timerName);
  if (timer.is_null()) {
    timer = TimeMonitor::getNewCounter(timerName);
  }
  return timer;
}

class FillSourceMatrixValues {
public:
  FillSourceMatrixValues(const crs_matrix_type::local_matrix_device_type& A,
                         const map_type::local_map_type& lclColMap)
    : A_(A), lclColMap_(lclColMap)
  {}
  KOKKOS_INLINE_FUNCTION void operator() (const LO lclRow) const {
    auto A_row = A_.row(lclRow);
    const double multiplier = 0.0; // (double(A_.numRows()) + 1.0) * double(lclRow);
    for(LO k = 0; k < A_row.length; ++k) {
      const LO lclCol = A_row.colidx(k);
      const GO gblCol = lclColMap_.getGlobalElement(lclCol);
      const double modVal = multiplier + (1.0 + double(gblCol));
      A_row.value(k) = modVal;
    }
  }
private:
  crs_matrix_type::local_matrix_device_type A_;
  map_type::local_map_type lclColMap_;
};

void
fillSourceMatrixValues(const crs_matrix_type& A)
{
  TEUCHOS_ASSERT(! A.isFillComplete() );
  auto A_lcl = A.getLocalMatrixDevice();

  TEUCHOS_ASSERT( A.hasColMap() );
  const auto lclColMap = A.getColMap()->getLocalMap();

  using execution_space =
    crs_matrix_type::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  Kokkos::parallel_for("Fill source's values",
                       range_type(0, A_lcl.numRows()),
                       FillSourceMatrixValues(A_lcl, lclColMap));
}

// FIXME (mfh 17 Nov 2019) We actually should have the source matrix
// have fewer entries per row than the target matrix, so that we don't
// hit any potential "all the column indices are the same"
// optimizations.

class TestTargetMatrixValues {
  using local_matrix_device_type = crs_matrix_type::local_matrix_device_type;
  using device_type = crs_matrix_type::device_type;
  using local_map_type = map_type::local_map_type;

public:
  TestTargetMatrixValues(const local_matrix_device_type& A,
                         const local_map_type& lclColMap,
                         const Kokkos::View<const int*, device_type>& lclRowsToTest,
                         const bool replaceCombineMode)
    : A_(A),
      lclColMap_(lclColMap),
      lclRowsToTest_(lclRowsToTest),
      replaceCombineMode_(replaceCombineMode)
  {}

  KOKKOS_INLINE_FUNCTION void init(int& success) const {
    success = 1;
  }

  KOKKOS_INLINE_FUNCTION void join(int& dst, int& src) const {
    dst = (dst == 1 && src == 1) ? 1 : 0;
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO lclRow, int& success) const {
    const bool modifiedRow = lclRowsToTest_[lclRow] != 0;

    auto A_row = A_.row(lclRow);
    const double multiplier = 0.0; // (double(A_.numRows()) + 1.0) * double(lclRow);
    for(LO k = 0; k < A_row.length; ++k) {
      const LO lclCol = A_row.colidx(k);
      const GO gblCol = lclColMap_.getGlobalElement(lclCol);
      const double modVal = multiplier + (1.0 + double(gblCol));
      const double modVal2 = replaceCombineMode_ ? modVal : modVal - 1.0;
      const double expVal = modifiedRow ? modVal2 : -1.0;

      const bool mySuccess = A_row.value(k) == expVal;
      success = mySuccess ? success : 0;
    }
  }
private:
  local_matrix_device_type A_;
  local_map_type lclColMap_;
  Kokkos::View<const int*, device_type> lclRowsToTest_;
  bool replaceCombineMode_;
};

int
testTargetLocalMatrixValues(const crs_matrix_type::local_matrix_device_type& A,
                            const map_type::local_map_type& lclColMap,
                            const Kokkos::View<const int*, typename crs_matrix_type::device_type>& lclRowsToTest,
                            const Tpetra::CombineMode combineMode,
                            const Teuchos::Comm<int>& comm,
                            const bool verbose)
{
  using std::endl;
  using execution_space =
    crs_matrix_type::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;

  const int myRank = comm.getRank();
  int lclSuccess = 1;

  const LO lclNumRows = A.numRows();
  // Kokkos::parallel_reduce("Test target's values",
  //   range_type(0, A.numRows()),
  //   TestTargetMatrixValues(A, lclColMap, lclRowsToTest, combineMode),
  //   Kokkos::Min<int, Kokkos::AnonymousSpace>(lclSuccess));
  const bool isReplace = combineMode == Tpetra::REPLACE;
  Kokkos::parallel_reduce("Test target's values",
    range_type(0, A.numRows()),
    TestTargetMatrixValues(A, lclColMap, lclRowsToTest, isReplace),
    lclSuccess);
  if(verbose) {
    std::ostringstream os_lcl;
    os_lcl << "Proc " << myRank << ": ";
    if(lclSuccess != 0) {
      os_lcl << "Local target matrix appears correct: "
        "lclSuccess=" << lclSuccess << endl;
    }
    else if(lclSuccess == 0) {
      os_lcl << "Local target matrix appears wrong.  "
             << "Let's retest sequentially on host." << endl;

      Kokkos::HostSpace hostMemSpace;
      using Kokkos::create_mirror_view;
      auto val = create_mirror_view(hostMemSpace, A.values);
      Kokkos::deep_copy(val, A.values);
      auto ind = create_mirror_view(hostMemSpace, A.graph.entries);
      Kokkos::deep_copy(ind, A.graph.entries);

      // create_mirror_view preserves const-ness of input View.
      using Kokkos::create_mirror;
      auto ptr = create_mirror(hostMemSpace, A.graph.row_map);
      Kokkos::deep_copy(ptr, A.graph.row_map);
      auto lclRowsToTest_h =
        create_mirror(hostMemSpace, lclRowsToTest);
      Kokkos::deep_copy(lclRowsToTest_h, lclRowsToTest);

      bool good = true;
      for(LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const bool modifiedRow = lclRowsToTest_h[lclRow];
        const LO numEnt = LO(ptr[lclRow+1] - ptr[lclRow]);

        const double multiplier = 0.0; // (double(A_.numRows()) + 1.0) * double(lclRow);
        for(LO k = 0; k < numEnt; ++k) {
          const LO lclCol = ind[k];
          const GO gblCol = lclColMap.getGlobalElement(lclCol);
          const double modVal = multiplier + (1.0 + double(gblCol));
          const double modVal2 = isReplace ? modVal : modVal - 1.0;
          const double expVal = modifiedRow ? modVal2 : -1.0;

          if(val[k] != expVal) {
            good = false;
            os_lcl << "Bad: lclRow=" << lclRow << ", lclCol="
                   << lclCol << ", val=" << val[k] << endl;
          }
        }
      }
      if(good) {
        os_lcl << "Local Matrix looks good here." << endl;
      }
    }

    if(myRank == 0) {
      std::cerr << "Local results:" << endl;
    }
    Tpetra::Details::gathervPrint(std::cerr, os_lcl.str(), comm);
    if(myRank == 0) {
      std::cerr << "Above were the local results." << endl;
    }
    comm.barrier();
  }
  return lclSuccess;
}

void
printLocalMatrix(std::ostream& out,
                 const crs_matrix_type& A)
{
  using std::endl;
  const auto& rowMap = *(A.getRowMap());
  auto comm = rowMap.getComm();
  const int myRank = comm->getRank();

  out << "Proc " << myRank << ":" << endl;

  auto A_lcl = A.getLocalMatrixHost();
  auto val = A_lcl.values;
  auto ind = A_lcl.graph.entries;
  auto ptr = A_lcl.graph.row_map;

  for(LO lclRow = 0; lclRow < A_lcl.numRows(); ++lclRow) {
    const GO gblRow = rowMap.getGlobalElement(lclRow);
    out << "(lclRow: " << lclRow << ", gblRow: " << gblRow << "): {";

    const ptrdiff_t beg = ptrdiff_t(ptr[lclRow]);
    const ptrdiff_t end = ptrdiff_t(ptr[lclRow+1]);
    for (ptrdiff_t k = beg; k < end; ++k) {
      out << "(" << ind[k] << "," << val[k] << ")";
      if (k + 1 < end) {
        out << ", ";
      }
    }
    out << "}" << endl;
  }
}

void
printMatrix(std::ostream& out,
            const crs_matrix_type& A,
            const std::string& label)
{
  using std::endl;

  std::ostringstream os_lcl;
  printLocalMatrix(os_lcl, A);
  std::ostringstream os_gbl;
  auto comm = A.getMap()->getComm();
  Tpetra::Details::gathervPrint(os_gbl, os_lcl.str(), *comm);
  if(comm->getRank() == 0) {
    out << label << ":" << endl << os_gbl.str() << endl;
  }
  // This is not strictly necessary, but it helps to let printing
  // finish before a test fails.
  comm->barrier();
}

bool
testTargetMatrixValues(const crs_matrix_type& A,
                       const Kokkos::View<const int*, typename crs_matrix_type::device_type>& lclRowsToTest,
                       const Tpetra::CombineMode combineMode)
{
  using std::endl;
  const bool verbose = Tpetra::Details::Behavior::verbose();
  auto comm = A.getMap()->getComm();
  const int myRank = comm->getRank();

  auto A_lcl = A.getLocalMatrixDevice();
  TEUCHOS_ASSERT( A.hasColMap() );
  auto lclColMap = A.getColMap()->getLocalMap();
  const int lclSuccess =
    testTargetLocalMatrixValues(A_lcl, lclColMap, lclRowsToTest,
                                combineMode, *comm, verbose);
  if(verbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": lclSuccess=" << lclSuccess << endl;
    std::cerr << os.str();
  }

  int gblSuccess = 1;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, lclSuccess,
                     Teuchos::outArg(gblSuccess));
  if(gblSuccess == 0) {
    std::ostringstream err;
    if(myRank == 0) {
      err << "Test failed: "
          << Tpetra::combineModeToString(combineMode)
          << " CombineMode" << endl;
    }
    return false;
  }
  else {
    return true;
  }
}

void benchmark(const RCP<const Teuchos::Comm<int>>& comm,
               const CmdLineArgs& args)
{
  using std::cerr;
  using std::endl;

  const int numProcs = comm->getSize();
  const int myRank = comm->getRank();
  const bool verbose = Tpetra::Details::Behavior::verbose();
  std::unique_ptr<std::string> prefix;
  if(verbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "benchmark: Make Maps" << endl;
    cerr << os.str();
  }

  Tpetra::Details::Behavior::disable_verbose_behavior();
  auto tgtRowMap = getTargetRowMap(comm, args);
  auto domainMap = getDomainMap(comm, args);
  auto rangeMap  = tgtRowMap;
  auto srcRowMap = getSourceRowMap(tgtRowMap, args);
  Tpetra::Details::Behavior::enable_verbose_behavior();

  RCP<Teuchos::FancyOStream> fancyOutPtr;
  if(verbose) {
    fancyOutPtr = Teuchos::getFancyOStream(Teuchos::rcpFromRef(cerr));
    cerr << "Target row Map:" << endl;
    tgtRowMap->describe(*fancyOutPtr, Teuchos::VERB_EXTREME);
    cerr << "Source row Map:" << endl;
    srcRowMap->describe(*fancyOutPtr, Teuchos::VERB_EXTREME);
  }

  const LO lclNumCols = LO(args.lclNumCols);
  const GO gblNumCols = GO(lclNumCols) * GO(numProcs);
  const LO maxNumColsToFill =
    std::ceil(args.densityFraction * double(gblNumCols));
  const LO numColsToFill = maxNumColsToFill > gblNumCols ?
    gblNumCols : maxNumColsToFill;

  if(verbose) {
    std::ostringstream os;
    os << *prefix << "Make Graphs: {lclNumCols: " << lclNumCols
       << ", numColsToFill: " << numColsToFill << "}" << endl;
    cerr << os.str();
  }

  Tpetra::Details::Behavior::disable_verbose_behavior();
  auto tgtGraph = rcp(new crs_graph_type(tgtRowMap, numColsToFill));
  auto srcGraph = rcp(new crs_graph_type(srcRowMap, numColsToFill));
  {
    Kokkos::fence(); // let device allocations & fills finish
    std::unique_ptr<GO[]> colGids(new GO[numColsToFill]);
    std::iota(colGids.get(), colGids.get()+numColsToFill,
              domainMap->getIndexBase());
    const LO tgtLclNumRows = LO(tgtRowMap->getLocalNumElements());
    for(LO lclRow = 0; lclRow < tgtLclNumRows; ++lclRow) {
      const GO gblRow = tgtRowMap->getGlobalElement(lclRow);
      tgtGraph->insertGlobalIndices(gblRow, numColsToFill,
                                    colGids.get());
    }
    const LO srcLclNumRows = LO(srcRowMap->getLocalNumElements());
    for(LO lclRow = 0; lclRow < srcLclNumRows; ++lclRow) {
      const GO gblRow = srcRowMap->getGlobalElement(lclRow);
      srcGraph->insertGlobalIndices(gblRow, numColsToFill,
                                    colGids.get());
    }
  }
  tgtGraph->fillComplete(domainMap, rangeMap);
  srcGraph->fillComplete(domainMap, rangeMap);

  crs_matrix_type tgtMatrix(tgtGraph);
  tgtMatrix.setAllToScalar(-1.0); // flag value
  tgtMatrix.fillComplete(domainMap, rangeMap);
  crs_matrix_type srcMatrix(srcGraph);

  fillSourceMatrixValues(srcMatrix);
  srcMatrix.fillComplete(domainMap, rangeMap);
  Tpetra::Details::Behavior::enable_verbose_behavior();

  if(verbose) {
    std::ostringstream os_lcl;
    printLocalMatrix(os_lcl, srcMatrix);
    std::ostringstream os_gbl;
    Tpetra::Details::gathervPrint(os_gbl, os_lcl.str(), *comm);
    if(myRank == 0) {
      cerr << "Source Matrix:" << endl << os_gbl.str() << endl;
    }
  }

  Tpetra::Details::Behavior::disable_verbose_behavior();
  export_type exporter(srcRowMap, tgtRowMap);
  Tpetra::Details::Behavior::enable_verbose_behavior();

  // Before benchmarking, actually test the two cases: CombineModes
  // REPLACE and ADD.

  {
    auto lclRowsToTest = getLclRowsToTest(exporter);
    tgtMatrix.resumeFill();
    {
      Tpetra::Details::Behavior::disable_verbose_behavior();
      const Tpetra::CombineMode combineMode = Tpetra::REPLACE;
      tgtMatrix.doExport(srcMatrix, exporter, combineMode);
      Tpetra::Details::Behavior::enable_verbose_behavior();
      if(verbose) {
        printMatrix(cerr, tgtMatrix, "Target Matrix (REPLACE)");
      }
      const bool ok =
        testTargetMatrixValues(tgtMatrix, lclRowsToTest, combineMode);
      if(! ok) {
        std::cerr << "Test FAILED" << std::endl;
      }
    }
    tgtMatrix.setAllToScalar(-1.0); // flag value
    {
      const Tpetra::CombineMode combineMode = Tpetra::ADD;
      tgtMatrix.doExport(srcMatrix, exporter, combineMode);
      if(verbose) {
        printMatrix(cerr, tgtMatrix, "Target Matrix (ADD)");
      }
      const bool ok =
        testTargetMatrixValues(tgtMatrix, lclRowsToTest, combineMode);
      if(! ok) {
        std::cerr << "Test FAILED" << std::endl;
      }
    }
    tgtMatrix.setAllToScalar(-1.0); // flag value
  }

  if(args.numTrials == 0) { // only test; don't benchmark
    return;
  }

  Tpetra::Details::Behavior::disable_verbose_behavior();
  for (auto combineMode : {Tpetra::ADD, Tpetra::REPLACE}) {
    auto timer = [&] {
      std::ostringstream os;
      os << "Tpetra::CrsMatrix::doExport: " << args.numTrials
         << " trial" << (args.numTrials != 1 ? "s" : "") << ", "
         << Tpetra::combineModeToString(combineMode) << " CombineMode";
      return getTimer(os.str());
    }();
    {
      TimeMonitor timeMon(*timer);
      for(int trial = 0; trial < args.numTrials; ++trial) {
        tgtMatrix.doExport(srcMatrix, exporter, Tpetra::ADD);
      }
    }
    tgtMatrix.setAllToScalar(-1.0); // flag value
  }
  Tpetra::Details::Behavior::enable_verbose_behavior();
}

void setCmdLineArgs(CmdLineArgs& args,
                    Teuchos::CommandLineProcessor& cmdp)
{
  cmdp.setOption ("densityFraction", &args.densityFraction,
                  "Fraction: Number of entries in each \"dense\" "
                  "row, divided by number of columns");
  cmdp.setOption ("lclNumRows", &args.lclNumRows,
                  "Number of row (and range) Map rows per process");
  cmdp.setOption ("lclNumCols", &args.lclNumCols,
                  "Number of domain Map indices per process; "
                  "-1 means \"same as lclNumRows\"");
  cmdp.setOption ("numTrials", &args.numTrials,
                  "Number of times to repeat each "
                  "operation in a timing loop");
  cmdp.setOption ("numOverlapRows", &args.numOverlapRows,
                  "Number of overlapping rows (that is, rows that "
                  "the target matrix needs to receive from a remote "
                  "process) per process.");
  cmdp.setOption ("contiguousMaps", "noncontiguousMaps",
                  &args.contiguousMaps,
                  "Whether the row, domain, and range Maps "
                  "of the target CrsMatrix are contiguous.");
}

enum class CmdLineArgsValidation {
  VALID,
  PRINTED_HELP,
  INVALID
};

CmdLineArgsValidation
getAndValidateCmdLineArgs(CmdLineArgs& args,
                          Teuchos::CommandLineProcessor& cmdp,
                          std::ostream& err,
                          int argc,
                          char* argv[])
{
  const CommandLineProcessor::EParseCommandLineReturn parseResult =
    cmdp.parse(argc, argv);
  if(parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    // The user specified --help at the command line to print help
    // with command-line arguments.
    return CmdLineArgsValidation::PRINTED_HELP;
  }
  else {
    using std::endl;
    bool valid = true;
    if(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL) {
      err << "Failed to parse command-line arguments." << endl;
      valid = false;
    }
    if(args.numTrials < 0) {
      err << "numTrials must be nonnegative." << endl;
      valid = false;
    }
    if(args.lclNumRows < 0) {
      err << "lclNumRows must be nonnegative." << endl;
      valid = false;
    }

    if(args.lclNumCols == -1) {
      args.lclNumCols = args.lclNumRows;
    }
    else if(args.lclNumCols < 0) {
      err << "lclNumCols must be nonnegative." << endl;
      valid = false;
    }

    if(args.numOverlapRows < 0) {
      err << "numOverlapRows must be nonnegative." << endl;
      valid = false;
    }
    if(args.numOverlapRows > args.lclNumRows) {
      err << "numOverlapRows must be <= lclNumRows." << endl;
      valid = false;
    }
    if(args.densityFraction < 0.0 || args.densityFraction > 1.0) {
      err << "densityFraction must be in [0,1]." << endl;
      valid = false;
    }
    return valid ? CmdLineArgsValidation::VALID :
      CmdLineArgsValidation::INVALID;
  }
}

void printCmdLineArgs(std::ostream& out,
                      const CmdLineArgs& args)
{
  using std::endl;
  out << endl << "---" << endl
      << "Command-line options:" << endl
      << "  densityFraction: " << args.densityFraction << endl
      << "  lclNumRows: " << args.lclNumRows << endl
      << "  lclNumCols: " << args.lclNumCols << endl
      << "  numTrials: " << args.numTrials << endl
      << "  numOverlapRows: " << args.numOverlapRows << endl
      << "  contiguousMaps: " << args.contiguousMaps << endl
      << endl;
}

int main(int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    auto tpetraComm = Tpetra::getDefaultComm();
    const int minNumProcs = 2;
    if(tpetraComm->getSize() < minNumProcs) {
      std::cerr << "This benchmark must be run with at least "
                << minNumProcs << " MPI processes." << std::endl;
      return EXIT_FAILURE;
    }

    CommandLineProcessor cmdp;
    CmdLineArgs args;
    setCmdLineArgs(args, cmdp);
    const auto validation =
      getAndValidateCmdLineArgs(args, cmdp, std::cerr, argc, argv);
    if(validation == CmdLineArgsValidation::PRINTED_HELP) {
      return EXIT_SUCCESS;
    }
    else if(validation == CmdLineArgsValidation::INVALID) {
      return EXIT_FAILURE;
    }
    else {
      if(tpetraComm->getRank() == 0) {
        printCmdLineArgs(std::cout, args);
      }
      benchmark(tpetraComm, args);
      TimeMonitor::report(tpetraComm.ptr(), std::cout);
    }
  }
  return EXIT_SUCCESS;
}
