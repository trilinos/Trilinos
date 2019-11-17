// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <numeric>
#include <cmath>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Time;
using Teuchos::TimeMonitor;

using GST = Tpetra::global_size_t;
using map_type = Tpetra::Map<>;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;
using crs_graph_type = Tpetra::CrsGraph<>;
using crs_matrix_type = Tpetra::CrsMatrix<double>;
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
  return std::move(tgtGids);
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
  FillSourceMatrixValues(const crs_matrix_type::local_matrix_type& A)
    : A_(A)
  {}
  KOKKOS_INLINE_FUNCTION void operator() (const LO lclRow) const {
    auto A_row = A_.row(lclRow);
    const double multiplier = (double(A_.numRows()) + 1.0) * double(lclRow);
    for(LO lclCol = 0; lclCol < A_row.length; ++lclCol) {
      A_row.value(lclCol) = multiplier + (1.0 + double(lclCol));
    }
  }
private:
  crs_matrix_type::local_matrix_type A_;
};

void
fillSourceMatrixValues(const crs_matrix_type& A)
{
  TEUCHOS_ASSERT(! A.isFillComplete() );
  auto A_lcl = A.getLocalMatrix();
  using execution_space =
    crs_matrix_type::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  Kokkos::parallel_for("Fill source's values",
                       range_type(0, A_lcl.numRows()),
                       FillSourceMatrixValues(A_lcl));
}

// FIXME (mfh 17 Nov 2019) We actually should have the source matrix
// have fewer entries per row than the target matrix, so that we don't
// hit any potential "all the column indices are the same"
// optimizations.

class TestTargetMatrixValues {
public:
  TestTargetMatrixValues(const crs_matrix_type::local_matrix_type& A,
                         const bool replaceCombineMode)
    : A_(A), replaceCombineMode_(replaceCombineMode)
  {}
  KOKKOS_INLINE_FUNCTION void
  operator() (const LO lclRow, int& success) const {
    auto A_row = A_.row(lclRow);
    const double multiplier = (double(A_.numRows()) + 1.0) * double(lclRow);
    for(LO lclCol = 0; lclCol < A_row.length; ++lclCol) {
      const double value = multiplier + (1.0 + double(lclCol));
      const bool mySuccess = replaceCombineMode_ ?
        value == A_row.value(lclCol) :
        value - 1.0 == A_row.value(lclCol);
      success = mySuccess ? success : 0;
    }
  }
private:
  crs_matrix_type::local_matrix_type A_;
  bool replaceCombineMode_;
};

int
testTargetLocalMatrixValues(const crs_matrix_type::local_matrix_type& A,
                            const Tpetra::CombineMode combineMode)
{
  const bool replaceCombineMode = combineMode == Tpetra::REPLACE;
  using execution_space =
    crs_matrix_type::device_type::execution_space;
  int lclSuccess = 0;
  Kokkos::parallel_reduce("Test target's values",
    Kokkos::RangePolicy<execution_space, LO>(0, A.numRows()),
    TestTargetMatrixValues(A, combineMode),
    Kokkos::Min<int, Kokkos::AnonymousSpace>(lclSuccess));
  return lclSuccess;
}

void
testTargetMatrixValues(const crs_matrix_type& A,
                       const Tpetra::CombineMode combineMode)
{
  const int lclSuccess =
    testTargetLocalMatrixValues(A.getLocalMatrix(), combineMode);
  int gblSuccess = 1;
  auto comm = A.getMap()->getComm();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, lclSuccess,
                     Teuchos::outArg(gblSuccess));
  TEUCHOS_TEST_FOR_EXCEPTION
    (gblSuccess != 1, std::logic_error, "Test failed: "
     << Tpetra::combineModeToString(combineMode) << " CombineMode");
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

  auto tgtRowMap = getTargetRowMap(comm, args);
  auto domainMap = getDomainMap(comm, args);
  auto rangeMap  = tgtRowMap;
  auto srcRowMap = getSourceRowMap(tgtRowMap, args);

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

  std::unique_ptr<GO[]> colGids(new GO[numColsToFill]);
  std::iota(colGids.get(), colGids.get()+numColsToFill,
            domainMap->getIndexBase());
  auto tgtGraph = rcp(new crs_graph_type(tgtRowMap, numColsToFill));
  auto srcGraph = rcp(new crs_graph_type(srcRowMap, numColsToFill));

  Kokkos::fence(); // let device allocations & fills finish
  for(LO lclRow = 0; lclRow < LO(tgtRowMap->getNodeNumElements()); ++lclRow) {
    const GO gblRow = tgtRowMap->getGlobalElement(lclRow);
    tgtGraph->insertGlobalIndices(gblRow, numColsToFill, colGids.get());
  }
  for(LO lclRow = 0; lclRow < LO(srcRowMap->getNodeNumElements()); ++lclRow) {
    const GO gblRow = srcRowMap->getGlobalElement(lclRow);
    srcGraph->insertGlobalIndices(gblRow, numColsToFill, colGids.get());
  }
  tgtGraph->fillComplete(domainMap, rangeMap);
  srcGraph->fillComplete(domainMap, rangeMap);

  crs_matrix_type tgtMatrix(tgtGraph);
  tgtMatrix.setAllToScalar(-1.0); // flag value
  tgtMatrix.fillComplete(domainMap, rangeMap);
  crs_matrix_type srcMatrix(srcGraph);

  fillSourceMatrixValues(srcMatrix);
  srcMatrix.fillComplete(domainMap, rangeMap);

  export_type exporter(srcRowMap, tgtRowMap);

  // Before benchmarking, actually test the two cases: CombineModes
  // REPLACE and ADD.

  tgtMatrix.resumeFill();
  {
    const Tpetra::CombineMode combineMode = Tpetra::REPLACE;
    tgtMatrix.doExport(srcMatrix, exporter, combineMode);
    testTargetMatrixValues(tgtMatrix, combineMode);
  }
  tgtMatrix.setAllToScalar(-1.0); // flag value
  {
    const Tpetra::CombineMode combineMode = Tpetra::ADD;
    tgtMatrix.doExport(srcMatrix, exporter, combineMode);
    testTargetMatrixValues(tgtMatrix, combineMode);
  }
  tgtMatrix.setAllToScalar(-1.0); // flag value

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
      << endl;
}

int main(int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    auto tpetraComm = Tpetra::getDefaultComm();
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
