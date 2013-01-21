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

#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_SingleLevelFactoryBase.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_AggOptions.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); // Note:

  GO nx,ny,nz;
  nx=500;
  ny=500;
  nz=100;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "scalingTest.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");
  int amgAsPrecond=1; clp.setOption("precond",&amgAsPrecond,"apply multigrid as preconditioner");
  int amgAsSolver=0; clp.setOption("fixPoint",&amgAsSolver,"apply multigrid as solver");
  bool printTimings=true; clp.setOption("timings","notimings",&printTimings,"print timings to screen");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (comm->getRank() == 0) {
    std::cout << "========================================================" << std::endl
              << xpetraParameters << matrixParameters;
  }

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));

  // read aggregation options from file
  Teuchos::FileInputSource fileSrc(xmlFileName);
  Teuchos::XMLObject fileXML = fileSrc.getObject();
  Teuchos::XMLParameterListReader listReader;
  Teuchos::ParameterList aggList = listReader.toParameterList(fileXML);
  //std::cout << "===========aggList start===========" << std::endl;
  //aggList.print(std::cout);
  //std::cout << "===========aggList end===========" << std::endl;

  // instantiate aggregate factory, set options from parameter list
  RCP<MueLu::SingleLevelFactoryBase> aggFact;
  if (aggList.name() == "UncoupledAggregationFactory") {

     RCP<UncoupledAggregationFactory> ucFact = rcp( new UncoupledAggregationFactory() );
     //ucFact->SetParameterList(aggList);
     //FIXME hack until UCAgg uses PL interface
     std::string ordering = aggList.get<std::string>("Ordering");
     MueLu::AggOptions::Ordering eordering;
     if (ordering=="Natural") eordering = MueLu::AggOptions::NATURAL;
     if (ordering=="Graph") eordering = MueLu::AggOptions::GRAPH;
     if (ordering=="Random") eordering = MueLu::AggOptions::RANDOM;
     ucFact->SetOrdering(eordering);
     ucFact->SetMaxNeighAlreadySelected(aggList.get<int>("MaxNeighAlreadySelected"));
     ucFact->SetMinNodesPerAggregate(aggList.get<int>("MinNodesPerAggregate"));
     aggFact = ucFact;

  } else if (aggList.name() == "CoupledAggregationFactory") {

     RCP<CoupledAggregationFactory> cFact = rcp( new CoupledAggregationFactory() );
     //cFact->SetParameterList(aggList);
     //FIXME hack until CoupledAgg uses PL interface
     //cFact->SetOrdering(aggList.get<std::string>("Ordering"));
     cFact->SetMaxNeighAlreadySelected(aggList.get<int>("MaxNeighAlreadySelected"));
     cFact->SetMinNodesPerAggregate(aggList.get<int>("MinNodesPerAggregate"));
     aggFact = cFact;

  } else {

    throw(MueLu::Exceptions::RuntimeError("List's name does not correspond to a known aggregation factory."));

  }

  //Teuchos::ParameterList tlist = aggFact->GetParameterList();
  //std::cout << "===========verify List start===========" << std::endl;
  //tlist.print(std::cout);
  //std::cout << "===========verify List end===========" << std::endl;

  // build matrix
  RCP<TimeMonitor> tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  RCP<const Map> map;
  RCP<MultiVector> coordinates;

  // Retrieve matrix parameters (they may have been changed on the command line), and pass them to Galeri.
  // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
  //                                 d1  d2  d3
  //                                 d4  d5  d6
  //                                 d7  d8  d9
  //                                 d10 d11 d12
  // A perfect distribution is only possible when the #processors is a perfect square.
  // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
  // size. For example, np=14 will give a 7-by-2 distribution.
  // If you don't want Galeri to do this, specify mx or my on the galeriList.
  Teuchos::ParameterList pl = matrixParameters.GetParameterList();
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", pl.get("nx",nx));
  galeriList.set("ny", pl.get("ny",ny));
  //galeriList.set("mx", comm->getSize());
  //galeriList.set("my", 1);

  if (matrixParameters.GetMatrixType() == "Laplace1D") {
    map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",map,matrixParameters.GetParameterList());
  }
  else if (matrixParameters.GetMatrixType() == "Laplace2D" || matrixParameters.GetMatrixType() == "Star2D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",map,matrixParameters.GetParameterList());
  }
  else if (matrixParameters.GetMatrixType() == "Laplace3D") {
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D",map,matrixParameters.GetParameterList());
    //map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList); //TODO when available in Galeri
    map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  }

  if (comm->getRank() == 0) {
    GO mx = galeriList.get("mx", -1);
    GO my = galeriList.get("my", -1);
    std::cout << "Processor subdomains in x direction: " << mx << std::endl
              << "Processor subdomains in y direction: " << my << std::endl
              << "========================================================" << std::endl;
  }

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
  RCP<Matrix> A = Pr->BuildMatrix();

  tm = Teuchos::null;

  Level level;
  RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  level.SetFactoryManager(factoryHandler);
  level.SetLevelID(0);
  level.Set("A", A);

  level.Request("Aggregates", aggFact.get());
  level.Request(*aggFact);

  level.setVerbLevel(Teuchos::VERB_NONE);
  aggFact->setVerbLevel(Teuchos::VERB_NONE);
  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("aggregation time")));
  aggFact->Build(level);
  tm = Teuchos::null;

  globalTimeMonitor = Teuchos::null;

  if (printTimings)
    TimeMonitor::summarize();

} //main
