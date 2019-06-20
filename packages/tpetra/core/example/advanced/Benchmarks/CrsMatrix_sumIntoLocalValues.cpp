/*
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <cstdlib>
#include <numeric>
#include <vector>

namespace { // (anonymous)

using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using TM = Teuchos::TimeMonitor;
using std::endl;
using crs_matrix_type = Tpetra::CrsMatrix<>;
using map_type = Tpetra::Map<>;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;

RCP<const Tpetra::Map<> >
createRowAndColMap (RCP<const Teuchos::Comm<int> > comm,
                    const LO lclNumInds)
{
  TM mon (*TM::getNewCounter ("createRowAndColMap"));

  const int numProcs = comm->getSize ();
  const GO gblNumInds = GO (numProcs) * GO (lclNumInds);
  const GO indexBase = 0;
  RCP<const map_type> rowAndColMap
    (new map_type (gblNumInds, lclNumInds, indexBase, comm));
  return rowAndColMap;
}

RCP<Tpetra::CrsMatrix<> >
createCrsMatrix (RCP<const map_type> rowAndColMap,
                 const size_t maxNumEntPerRow)
{
  TM mon (*TM::getNewCounter ("CrsMatrix constructor"));

  return rcp (new crs_matrix_type (rowAndColMap, rowAndColMap,
                                   maxNumEntPerRow,
                                   Tpetra::StaticProfile));
}

void
populateCrsMatrix (crs_matrix_type& A,
                   const Tpetra::CrsMatrix<>::local_ordinal_type numToInsert,
                   const Tpetra::CrsMatrix<>::local_ordinal_type lclColInds[],
                   const double vals[])
{
  TM mon (*TM::getNewCounter ("insertLocalValues loop"));

  const LO lclNumRows (A.getNodeNumRows ());
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    A.insertLocalValues (lclRow, numToInsert, vals, lclColInds);
  }
}

void
doSumIntoLocalValues (const std::string& label,
                      crs_matrix_type& A,
                      const LO numToInsert,
                      const LO lclColInds[],
                      const double vals[],
                      const int numTrials)
{
  TM mon (*TM::getNewCounter (label));

  const LO lclNumRows (A.getNodeNumRows ());

  for (int trial = 0; trial < numTrials; ++trial) {
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      (void) A.sumIntoLocalValues (lclRow, numToInsert, vals, lclColInds);
    }
  }
}

struct CmdLineArgs {
  int lclNumInds = 10000;
  int maxNumEntPerRow = 28;
  int numToInsert = 27;
  int numTrials = 10;
};

void
benchmarkCrsMatrixSumIntoLocalValues (const CmdLineArgs args)
{
  TM mon (*TM::getNewCounter ("Whole benchmark"));

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  auto rowAndColMap = createRowAndColMap (comm, LO (args.lclNumInds));
  auto A = createCrsMatrix (rowAndColMap, size_t (args.maxNumEntPerRow));

  std::vector<LO> lclColInds (size_t (args.numToInsert));
  std::vector<double> vals (size_t (args.numToInsert), 0.0);
  // Skip 0, so that search isn't completely trivial.
  std::iota (lclColInds.begin (), lclColInds.end (), LO (1));

  populateCrsMatrix (*A, LO (args.numToInsert),
                     lclColInds.data (), vals.data ());
  doSumIntoLocalValues ("sumIntoLocalValues: before first fillComplete",
                        *A, LO (args.numToInsert), lclColInds.data (),
                        vals.data (), args.numTrials);
  {
    TM mon2 (*TM::getNewCounter ("CrsMatrix fillComplete"));
    A->fillComplete ();
  }
  doSumIntoLocalValues ("sumIntoLocalValues: after first fillComplete",
                        *A, LO (args.numToInsert), lclColInds.data (),
                        vals.data (), args.numTrials);
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using Teuchos::CommandLineProcessor;
  using std::cerr;
  using std::cout;
  using std::endl;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  CommandLineProcessor cmdp;
  CmdLineArgs cmdLineArgs;
  cmdp.setOption ("lclNumInds", &cmdLineArgs.lclNumInds,
                  "Number of rows on each (MPI) process");
  cmdp.setOption ("maxNumEntPerRow", &cmdLineArgs.maxNumEntPerRow,
                  "Max amount of space to preallocate per row "
                  "(must be at least as big as numToInsert)");
  cmdp.setOption ("numToInsert", &cmdLineArgs.numToInsert,
                  "Number of column indices to insert per row");
  cmdp.setOption ("numTrials", &cmdLineArgs.numTrials,
                  "Number of times to repeat each operation in a "
                  "timing loop");
  const CommandLineProcessor::EParseCommandLineReturn parseResult =
    cmdp.parse (argc, argv);
  if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    // The user specified --help at the command line to print help
    // with command-line arguments.  We printed help already, so quit
    // with a happy return code.
    return EXIT_SUCCESS;
  }
  else if (parseResult != CommandLineProcessor::PARSE_SUCCESSFUL) {
    cerr << "Failed to parse command-line arguments." << endl;
    return EXIT_FAILURE;
  }
  cout << "lclNumInds: " << cmdLineArgs.lclNumInds << endl
       << "maxNumEntPerRow: " << cmdLineArgs.maxNumEntPerRow << endl
       << "numToInsert: " << cmdLineArgs.numToInsert << endl
       << "numTrials: " << cmdLineArgs.numTrials << endl;
  benchmarkCrsMatrixSumIntoLocalValues (cmdLineArgs);

  auto comm = Tpetra::getDefaultComm ();
  TM::report (comm.ptr (), cout);

  return EXIT_SUCCESS;
}
