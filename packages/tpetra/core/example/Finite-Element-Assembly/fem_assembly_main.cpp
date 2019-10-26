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
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

#include "fem_assembly_commandLineOpts.hpp"
#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_InsertGlobalIndices_FE_SP.hpp"
#include "fem_assembly_InsertGlobalIndices_DP.hpp"
#include "fem_assembly_LocalElementLoop_DP.hpp"
#include "fem_assembly_TotalElementLoop_DP.hpp"
#include "fem_assembly_TotalElementLoop_SP.hpp"


using namespace TpetraExamples;


int main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::StackedTimer;

  int status = EXIT_SUCCESS;

  // MPI boilerplate
  Tpetra::initialize(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out = *pOut;

  // Read command-line options into the 'opts' struct.
  struct CmdLineOpts opts;

  try
  {
    status = readCmdLineOpts(out, opts, argc, argv);
  }
  catch(...)
  {
    status = EXIT_FAILURE;
  }

  if(EXIT_SUCCESS != status)
  {
    Tpetra::finalize();
    return status;
  }

  RCP<StackedTimer> timer = Teuchos::null;
  if(opts.timing)
  {
    timer = rcp(new StackedTimer("X) Global", false));
    TimeMonitor::setStackedTimer(timer);
  }

  // Entry point
  if(opts.useStaticProfile)
  {
    if(opts.execInsertGlobalIndicesFE && executeInsertGlobalIndicesFESP(comm, opts))
       status = EXIT_FAILURE;
    if(opts.execTotalElementLoop && executeTotalElementLoopSP(comm, opts))
      status = EXIT_FAILURE;
  }
  else
  {
    if(opts.execInsertGlobalIndices && executeInsertGlobalIndicesDP(comm, opts))
       status = EXIT_FAILURE;
    if(opts.execLocalElementLoop    && executeLocalElementLoopDP(comm, opts))
       status = EXIT_FAILURE;
    if(opts.execTotalElementLoop    && executeTotalElementLoopDP(comm, opts))
       status = EXIT_FAILURE;
  }

  if(opts.timing)
  {
    //note: base timer was already stopped by executeInsertGlobalIndices...()
    StackedTimer::OutputOptions timeReportOpts;
    timeReportOpts.print_warnings = false;
    timer->report(std::cout, comm, timeReportOpts);
  }

  // This tells the Trilinos test framework that the test passed.
  if(EXIT_SUCCESS == comm->getRank()) out << "End Result: TEST PASSED" << endl;
  else                                out << "End Result: TEST FAILED" << endl;

  // Finalize
  Tpetra::finalize();

  return status;
}  // END main()

