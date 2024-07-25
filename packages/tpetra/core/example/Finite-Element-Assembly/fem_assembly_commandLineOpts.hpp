// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP

#include "Teuchos_CommandLineProcessor.hpp"

namespace TpetraExamples {

// Options to read in from the command line
struct CmdLineOpts
{
  // numElementsX - Number of Elements to generate in the X dimension.
  size_t numElementsX;
  // numElementsY - Number of Elements to generate in the Y dimension.
  size_t numElementsY;
  // verbose - Verbose mode, print out information messages
  bool verbose;
  // timing - Print out timing information at the end of the run.
  bool timing;
  // saveMM - Save the crsMatrix into a MatrixMarket file.
  bool saveMM;
  // repetitions - how many times to execute the kernel for testing
  size_t repetitions;
  // Use Kokkos assembly for matrix
  bool useKokkosAssembly;
  // Number of doubles per element (for simulated state transfer)
  int numStateDoublesPerElement;
  // execInsertGlobalIndicesFE - execute the FE Insert Global Indices kernel
  bool execInsertGlobalIndicesFE;
  // execTotalElementLoop - execute the Total Element Loop kernel
  bool execTotalElementLoop;
};

// Use a utility from the Teuchos package of Trilinos to set up
// command-line options for reading, and set default values of
// command-line options.  clp is an output argument containing the
// set-up options.  It retains pointers to fields in 'opts'.
// Reading the command-line options will update those fields in
// place.
void
setCmdLineOpts (struct CmdLineOpts& opts,
                Teuchos::CommandLineProcessor& clp)
{
  // Set default values of command-line options.
  opts.numElementsX = 3;
  opts.numElementsY = 3;
  opts.verbose = false;
  opts.timing = true;
  opts.saveMM = false;
  opts.repetitions  = 1;
  opts.useKokkosAssembly = false;
  opts.numStateDoublesPerElement = 4;
  opts.execInsertGlobalIndicesFE = false;
  opts.execTotalElementLoop = false;

  clp.setOption("num-elements-x", &(opts.numElementsX), "Number of elements to generate in the X-directon of the 2D grid.");
  clp.setOption("num-elements-y", &(opts.numElementsY), "Number of elements to generate in the Y-direction of the 2D grid.");

  clp.setOption("verbose", "without-verbose", &(opts.verbose), "Execute example with high verbosity.");
  clp.setOption("timing",  "without-timing",  &(opts.timing),  "Print out timing information at the end of execution.");
  clp.setOption("save-mm", "without-save-mm", &(opts.saveMM),  "Save the generated CrsMatrix into a Matrix Market file.");
  clp.setOption("repetitions", &(opts.repetitions), "Number of times to repeat the kernel.");
  clp.setOption("kokkos", "no-kokkos", &(opts.useKokkosAssembly), "Use Kokkos assembly.");
  clp.setOption("state-per-element",&(opts.numStateDoublesPerElement),"Number of doubles per element to store element state");
  clp.setOption("with-insert-global-indices-fe", "without-insert-global-indices-fe", &(opts.execInsertGlobalIndicesFE),
                "Execute the Insert FECrsMatrix Global Indices FEM Assembly kernel.");
  clp.setOption("with-total-element-loop", "without-total-element-loop", &(opts.execTotalElementLoop),
                "Execute the Total Element Loop FEM Assembly kernel.");
}

// Check the command-line options that were read in by
// parseCmdLineOpts.  Return 0 if all correct, else return nonzero,
// using the LAPACK error reporting convention of the negative of
// the argument in its original order (starting with 1) as the error
// code.  Print informative error messages to the given output
// stream \c out.
int checkCmdLineOpts(std::ostream& out, const struct CmdLineOpts& opts)
{
  int err = 0;
  // Currently we only have StaticProfile for TotalElementLoop
  if(1 != (opts.execInsertGlobalIndicesFE + opts.execTotalElementLoop))
  {
    out << std::endl
      << "  The FE assembly has two main modes, insert-global-indices-fe and total-element-loop." << std::endl
      << "  Exactly one of them should be enabled in a given run." << std::endl
      << "  --with-insert-global-indices-fe : Execute the FE Insert Global Indices example." << std::endl
      << "  --with-total-element-loop       : Execute the Total Element Loop example." << std::endl
      << std::endl;
    err = -1;
  }
  return err;
}

// Actually read the command-line options from the command line,
// using the argc and argv arguments to main().  Use the clp output
// argument of setCmdLineOpts.  The options actually get written to
// the same CmdLineOpts struct instance that was passed into
// setCmdLineOpts above.
//
// Return 0 if successful, 1 if help printed due to the user
// specifying the "--help" option (indicates that the application
// shouldn't run), and -1 on error.
int parseCmdLineOpts(Teuchos::CommandLineProcessor& clp, int argc, char* argv[])
{
  auto result = clp.parse(argc, argv);

  switch(result)
  {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return 1;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return -1;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      return 0;
    default:
      return -1;
  }
}

int readCmdLineOpts(std::ostream& out, struct CmdLineOpts& opts, int argc, char* argv[])
{
  using std::endl;

  {
    Teuchos::CommandLineProcessor clp;
    setCmdLineOpts(opts, clp);
    int result = parseCmdLineOpts(clp, argc, argv);
    // help printed
    if(1 == result)
    {
      return EXIT_SUCCESS;
    }
    // parse error
    else if(-1 == result)
    {
      return EXIT_FAILURE;
    }
    result = checkCmdLineOpts(out, opts);
    if(0 != result)
    {
      return EXIT_FAILURE;
    }
  }

  if (opts.verbose) {
    out << "Command-line options:" << endl;

    Teuchos::OSTab tab1(out); // push one tab in this scope
    out << "numElementsX : " << opts.numElementsX     << endl
        << "numElementsY : " << opts.numElementsY     << endl
        << "verbose      : " << opts.verbose          << endl
        << "timing       : " << opts.timing           << endl
        << "saveMM       : " << opts.saveMM           << endl
        << "repetitions  : " << opts.repetitions      << endl
        << endl
        << "execInsertGlobalIndicesFE : " << opts.execInsertGlobalIndicesFE << endl
        << "execTotalElementLoop      : " << opts.execTotalElementLoop    << endl
        << endl;
  }

  return EXIT_SUCCESS;
}

} // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP

