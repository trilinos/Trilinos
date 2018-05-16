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
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP

#include <Teuchos_CommandLineProcessor.hpp>


namespace TpetraExamples
{


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
  // StaticProfile
  bool useStaticProfile;
  // execInsertGlobalIndicesDP - execute the Insert Global Indices kernel
  bool execInsertGlobalIndices;
  // execLocalElementLoopDP - execute the Local Element Loop kernel
  bool execLocalElementLoop;
  // execTotalElementLoopDP - execute the Total Element Loop kernel
  bool execTotalElementLoop;
  // repetitions - how many times to execute the kernel for testing
  size_t repetitions;
};



// Use a utility from the Teuchos package of Trilinos to set up
// command-line options for reading, and set default values of
// command-line options.  clp is an output argument containing the
// set-up options.  It retains pointers to fields in 'opts'.
// Reading the command-line options will update those fields in
// place.
void setCmdLineOpts(struct CmdLineOpts& opts, Teuchos::CommandLineProcessor& clp)
{
  // Set default values of command-line options.
  opts.numElementsX = 3;
  opts.numElementsY = 3;
  opts.verbose = false;
  opts.timing = true;
  opts.saveMM = false;
  opts.useStaticProfile  = true;
  opts.execInsertGlobalIndices = false;
  opts.execLocalElementLoop    = false;
  opts.execTotalElementLoop    = false;
  opts.repetitions  = 1;

  clp.setOption("num-elements-x", &(opts.numElementsX), "Number of elements to generate in the X-directon of the 2D grid.");
  clp.setOption("num-elements-y", &(opts.numElementsY), "Number of elements to generate in the Y-direction of the 2D grid.");

  clp.setOption("verbose", "without-verbose", &(opts.verbose), "Execute example with high verbosity.");
  clp.setOption("timing",  "without-timing",  &(opts.timing),  "Print out timing information at the end of execution.");
  clp.setOption("save-mm", "without-save-mm", &(opts.saveMM),  "Save the generated CrsMatrix into a Matrix Market file.");

  clp.setOption("with-StaticProfile", "with-DynamicProfile", &(opts.useStaticProfile), "Use StaticProfile or DynamicProfile");

  clp.setOption("with-insert-global-indices", "without-insert-global-indices", &(opts.execInsertGlobalIndices),
                "Execute the Insert Global Indices FEM Assembly kernel.");
  clp.setOption("with-local-element-loop",    "without-local-element-loop",    &(opts.execLocalElementLoop),
                "Execute the Local Element Loop FEM Assembly kernel.");
  clp.setOption("with-total-element-loop",    "without-total-element-loop",    &(opts.execTotalElementLoop),
                "Execute the Total Element Loop FEM Assembly kernel.");
  clp.setOption("repetitions", &(opts.repetitions), "Number of times to repeat the kernel.");
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



// Check the command-line options that were read in by
// parseCmdLineOpts.  Return 0 if all correct, else return nonzero,
// using the LAPACK error reporting convention of the negative of
// the argument in its original order (starting with 1) as the error
// code.  Print informative error messages to the given output
// stream \c out.
int checkCmdLineOpts(std::ostream& out, const struct CmdLineOpts& opts)
{
  int err = 0;

  if( 1 != (opts.execInsertGlobalIndices + opts.execLocalElementLoop + opts.execTotalElementLoop))
  {
    out << std::endl
        << "Please select one algorithm to run.  Options are:" << std::endl
        << "  --with-insert-global-indices  :  Execute the Insert Global Indices example." << std::endl
        << "  --with-local-element-loop     :  Execute the Local Element Loop example." << std::endl
        << "  --with-total-element-loop     :  Execute the Total Element Loop example." << std::endl
        << std::endl;
    err = -1;
  }
  else
  {
    // Currently we only have StaticProfile for TotalElementLoop
    if(opts.useStaticProfile && !(opts.execTotalElementLoop))
    {
      out << std::endl
          << "StaticProfile is currently only implemented with Total Element Loop, please use:" << std::endl
          << "  --with-total-element-loop :  Execute the Total Element Loop example." << std::endl
          << std::endl;
      err = 1;
    }
  }
  return err;
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

  if(opts.verbose)
  {
    out << "Command-line options:" << endl;
    {
      Teuchos::OSTab tab1(out); // push one tab in this scope
      out << "numElementsX : " << opts.numElementsX     << endl
          << "numElementsY : " << opts.numElementsY     << endl
          << "verbose      : " << opts.verbose          << endl
          << "timing       : " << opts.timing           << endl
          << "saveMM       : " << opts.saveMM           << endl
          << "staticProfile: " << opts.useStaticProfile << endl
          << "repetitions  : " << opts.repetitions      << endl
          << endl
          << "execInsertGlobalIndices: " << opts.execInsertGlobalIndices << endl
          << "execLocalElementLoop   : " << opts.execLocalElementLoop    << endl
          << "execTotalElementLoop   : " << opts.execTotalElementLoop    << endl
          << endl;
    }
  }

  return EXIT_SUCCESS;
}



} // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_COMMANDLINEOPTS_HPP

