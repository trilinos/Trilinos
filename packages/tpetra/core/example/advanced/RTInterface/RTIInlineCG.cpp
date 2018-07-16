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

#include "RTIInlineCG.hpp"
#include "Tpetra_Core.hpp"

int main (int argc, char *argv[])
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ParameterList;

  //
  // Get the communicator
  //
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();

  //
  // Get example parameters from command-line processor
  //
  bool verbose = (myRank==0);
  std::string xmlfile;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // instantiate a driver on the scalar stack
  //
  CGDriver<double> driver;
  // hand output stream to driver
  if (verbose) {
    driver.out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  } else {
    RCP<Teuchos::oblackholestream> blackhole =
      rcp (new Teuchos::oblackholestream ());
    driver.out = Teuchos::getFancyOStream (blackhole);
  }

  //
  // get the solver parameters
  //
  RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                                  \n"
    "   <Parameter name='tolerance' value='5e-15' type='double'/>      \n"
    "   <Parameter name='verbose' value='1' type='int'/>               \n"
    "   <Parameter name='kappa' value='100' type='double'/>            \n"
    "   <Parameter name='size'  value='100' type='int'/>               \n"
    "   <Parameter name='numIters' value='100' type='int'/>            \n"
    " </ParameterList>                                                 \n"
  );
  Teuchos::updateParametersFromXmlString (xmlString,params.ptr ());
  if (xmlfile != "") {
    Teuchos::updateParametersFromXmlFile (xmlfile, params.ptr ());
  }
  // hand solver parameters to driver
  driver.params = params;

  //
  // run the driver
  //
  driver.run (comm);

  //
  // Print result
  if (driver.testPassed) {
    *driver.out << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
