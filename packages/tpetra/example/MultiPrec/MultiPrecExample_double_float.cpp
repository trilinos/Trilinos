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

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_DefaultPlatform.hpp>

// no CXX11 on the GPU
#undef HAVE_KOKKOSCLASSIC_THRUST
#include <Tpetra_HybridPlatform.hpp>
#include <TpetraExt_TypeStack.hpp>

#include <iostream>

#include "MultiPrecDriver.hpp"

/** \file MultiPrecExample_double_float.cpp
    \brief An example of a multi-precision algorithm, using a flexible preconditioned CG with recursive precision.
 */

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ParameterList;

  // 
  // Get the communicator
  //
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myImageID = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  bool unfused = false;
  std::string matfile;
  std::string xmlfile;
  std::string machineFile;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("matrix-file",&matfile,"Filename for matrix");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  cmdp.setOption("machine-file",&machineFile,"Filename for XML machine description file.");
  cmdp.setOption("unfused","no-unfused",&unfused,"Test unfused iteration.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // read machine file and initialize platform
  // 
  RCP<Teuchos::ParameterList> machinePL = Teuchos::parameterList();
  std::string defaultMachine(
    " <ParameterList>                                                               "
    "   <ParameterList name='%1=0'>                                                 "
    "     <Parameter name='NodeType'     type='string' value='Kokkos::SerialNode'/> "
    "   </ParameterList>                                                            "
    " </ParameterList>                                                              "
  );
  Teuchos::updateParametersFromXmlString(defaultMachine,machinePL.ptr());
  if (machineFile != "") Teuchos::updateParametersFromXmlFile(machineFile,machinePL.ptr());

  // 
  // create the platform object
  // 
  Tpetra::HybridPlatform platform(comm,*machinePL);

  // 
  // Define the type stack
  // 
  TPETRAEXT_TYPESTACK2(MPStack, double, float)

  //
  // instantiate a driver on the scalar stack
  //
  MultiPrecDriver<MPStack> driver;
  // hand output stream to driver
  if (verbose) driver.out = Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false));
  else         driver.out = Teuchos::getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));
  // hand matrix file to driver
  driver.matrixFile = matfile;
  // other params
  driver.unfusedTest = unfused;

  //
  // get the solver parameters
  // 
  RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                                    \n"
    "   <Parameter name='tolerance' value='1e-15' type='double'/>        \n"
    "   <Parameter name='verbose' value='2' type='int'/>                 \n"
    "   <ParameterList name='child'>                                     \n"
    "     <Parameter name='tolerance' value='1e-7' type='double'/>       \n"
    "     <Parameter name='verbose' value='1' type='int'/>               \n"
    "     <Parameter name='Extract Diagonal' value='true' type='bool'/>  \n" 
    "   </ParameterList>                                                 \n"
    " </ParameterList>                                                   \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,params.ptr());
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,params.ptr());
  // hand solver parameters to driver
  driver.params = params;

  // 
  // run the driver
  // 
  platform.runUserCode(driver);

  //
  // Print result
  if (driver.testPassed) {
    *driver.out << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}

/** \example MultiPrecExample_double_float.cpp 
    Demonstrate using Tpetra::RTI and a multi-precision flexible preconditioned CG, Tpetra::TypeStack and related utilities.
  */
