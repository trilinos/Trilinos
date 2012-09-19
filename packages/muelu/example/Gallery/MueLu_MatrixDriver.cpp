// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp> 

/*
  This driver simply generates a Tpetra matrix, prints it to screen, and exits.

  THIS EXAMPLE DOES NOT USE XPETRA.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char** argv) 
{
  using Teuchos::RCP;

  typedef int    LO; // LocalOrdinal
  typedef int    GO; // GlobalOrdinal
  typedef double SC; // Scalar

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp);   // manage parameters of the test case
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  std::cout << matrixParameters;

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  RCP<const Tpetra::Map<LO,GO> > map = rcp( new Tpetra::Map<LO,GO>(matrixParameters.GetNumGlobalElements(), 0, comm) );
  RCP<Tpetra::CrsMatrix<SC,LO,GO> > A = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Tpetra::Map<LO,GO>, Tpetra::CrsMatrix<SC,LO,GO> >(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if (comm->getRank() == 0)
    std::cout << "\n================ MAP =====================================================\n" << std::endl;
  map->describe(*out, Teuchos::VERB_EXTREME);
  comm->barrier();
  sleep(1);

  if (comm->getRank() == 0)
    std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
  A->describe(*out, Teuchos::VERB_EXTREME);

  return EXIT_SUCCESS;
} 
