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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include "MueLu_SemiCoarsenPFactory.hpp" // for semi-coarsening constants

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

// These files must be included last
#include <MueLu_UseDefaultTypes.hpp>
#include <unistd.h>
/**********************************************************************************/

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  bool success = true;
  bool verbose = true;
  try {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    /**********************************************************************************/
    /* SET TEST PARAMETERS                                                            */
    /**********************************************************************************/
    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);

    // Default is Laplace1D with nx = 8748.
    // It's a nice size for 1D and perfect aggregation. (6561=3^8)
    //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

    // custom parameters
    int pauseForDebugger=0;
    //std::string aggOrdering = "natural";
    int minPerAgg=2; //was 3 in simple
    int maxNbrAlreadySelected=0;
    int printTimings=0;
    std::string xmlFile="parameters.xml";

    //clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,graph)");
    clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
    clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
    clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
    clp.setOption("timings",&printTimings,"print timings to screen");
    clp.setOption("xmlFile",&xmlFile,"file name containing MueLu multigrid parameters in XML format");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor = Teuchos::rcp (new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Timings: Global Time")));

    if (pauseForDebugger) {
      Utils::PauseForDebugger();
    }

    matrixParameters.check();
    xpetraParameters.check();
    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (comm->getRank() == 0) {
      std::cout << xpetraParameters << matrixParameters;
    }

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    Teuchos::RCP<const Map> map;
    Teuchos::RCP<Matrix> A;

    {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("Timings: Matrix Build"));

      map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
      Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap
      A = Pr->BuildMatrix();

    }
    /**********************************************************************************/
    /*                                                                                */
    /**********************************************************************************/
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

    // create parameter list interpreter
    Teuchos::RCP<HierarchyManager> mueluFactory = Teuchos::rcp(new ParameterListInterpreter(paramList));

    Teuchos::RCP<Hierarchy> H = mueluFactory->CreateHierarchy();

    H->GetLevel(0)->Set< Teuchos::RCP<Matrix> >("A", A);

    Teuchos::RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getRowMap(), 1);
    nullspace->putScalar(1.0);
    H->GetLevel(0)->Set("Nullspace", nullspace);

    // set minimal information about number of layers for semicoarsening...
    // This information can also be provided as a user parameter in the xml file using the
    // parameter: "semicoarsen: num layers"
    H->GetLevel(0)->Set("NumZLayers",matrixParameters.GetParameterList().get<GO>("nz"));


    mueluFactory->SetupHierarchy(*H);

    for (int l=0; l<H->GetNumLevels(); l++) {
      Teuchos::RCP<MueLu::Level> level = H->GetLevel(l);
      if(level->IsAvailable("A", MueLu::NoFactory::get()) == false) { success = false; }
      if(level->IsAvailable("P", MueLu::NoFactory::get()) == false && l>0) { success = false; }
      if(level->IsAvailable("R", MueLu::NoFactory::get()) == false && l>0) { success = false; }
      if(level->IsAvailable("PreSmoother",  MueLu::NoFactory::get()) == false) { success = false; }
      if(level->IsAvailable("PostSmoother", MueLu::NoFactory::get()) == false && l<H->GetNumLevels()-1) { success = false; }
      if(level->IsAvailable("NumZLayers",   MueLu::NoFactory::get()) == true) {  success = false; }
      //H->GetLevel(l)->print(std::cout, MueLu::Debug);
    }
    ///////////////////////////////////////////////////////////

    // Timer final summaries
    globalTimeMonitor = Teuchos::null; // stop this timer before summary

    if (printTimings)
      Teuchos::TimeMonitor::summarize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
