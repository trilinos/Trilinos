// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Tpetra_Details_KokkosTeuchosTimerInjection.hpp"
#include "fem_assembly_commandLineOpts.hpp"
#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_InsertGlobalIndices_FE.hpp"
#include "fem_assembly_TotalElementLoop.hpp"


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

  // Force timing of the Kokkos::deep_copy calls
  Tpetra::Details::AddKokkosDeepCopyToTimeMonitor(true);

  // Entry point
  if(opts.execInsertGlobalIndicesFE && executeInsertGlobalIndicesFESP(comm, opts))
     status = EXIT_FAILURE;
  if(opts.execTotalElementLoop && executeTotalElementLoopSP(comm, opts))
    status = EXIT_FAILURE;

  if(opts.timing)
  {
    //note: base timer was already stopped by executeInsertGlobalIndices...()
    StackedTimer::OutputOptions timeReportOpts;
    timeReportOpts.print_warnings = false;
    timer->report(std::cout, comm, timeReportOpts);
    auto xmlOut = timer->reportWatchrXML("Tpetra FE Assembly " + std::to_string(comm->getSize()) + " ranks", comm);
    if(xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
  }

  // This tells the Trilinos test framework that the test passed.
  if(EXIT_SUCCESS == comm->getRank()) out << "End Result: TEST PASSED" << endl;
  else                                out << "End Result: TEST FAILED" << endl;

  // Finalize
  Tpetra::finalize();

  return status;
}  // END main()

