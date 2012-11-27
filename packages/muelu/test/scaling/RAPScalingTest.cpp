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

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include "MueLu.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TransPFactory.hpp"

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>


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
  nx=350;
  ny=350;
  nz=350;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  int  numRAPs=100; clp.setOption("nraps",&numRAPs,"number of RAPS to perform");
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

  //
  // Construct the problem
  //

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
  RCP<TimeMonitor> tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("RAPScalingTest: 1 - Matrix creation")));

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
  galeriList.set("mx", comm->getSize()); //mx and my are the number of domains
  galeriList.set("my", 1);

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
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
  }

//FIXME
  if (comm->getRank() == 0) {
    GO mx = galeriList.get("mx", -1);
    GO my = galeriList.get("my", -1);
    std::cout << "Processor subdomains in x direction: " << mx << std::endl
              << "Processor subdomains in y direction: " << my << std::endl
              << "========================================================" << std::endl;
  }

  RCP<Matrix> A   = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("RAPScalingTest: 2 - Setup")));

  Level fineLevel, coarseLevel;
  MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy( fineLevel, coarseLevel);

  fineLevel.Set("A",A);

  TentativePFactory tentpFactory;
  SaPFactory sapFactory(rcpFromRef(tentpFactory));
  TransPFactory transPFactory;
  transPFactory.SetFactory("P", rcpFromRef(sapFactory));
  coarseLevel.Request("P", &sapFactory);
  coarseLevel.Request("R", &transPFactory);

  sapFactory.Build(fineLevel, coarseLevel);
  transPFactory.Build(fineLevel,coarseLevel);
  RAPFactory rap;
  rap.SetFactory("P", rcpFromRef(sapFactory));
  rap.SetFactory("R", rcpFromRef(transPFactory));

  for (int i=0; i<numRAPs; ++i) {
    coarseLevel.Request("A",&rap);
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("RAPScalingTest: 3 - RAP kernel")));
    RCP<Matrix> Ac = coarseLevel.Get< RCP<Matrix> >("A", &rap);
    tm = Teuchos::null;
    coarseLevel.Release("A",&rap);
  }

  tm = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (printTimings) {
    Teuchos::TableFormat &format = TimeMonitor::format();
    format.setPrecision(25);
    TimeMonitor::summarize();
  }

} //main
