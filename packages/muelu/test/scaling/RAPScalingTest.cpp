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

#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include "MueLu.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"

#include <MueLu_UseDefaultTypes.hpp>

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false);

  GO nx = 50, ny = 50, nz = 50;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName       = "rapTest.xml"; clp.setOption("xml",                   &xmlFileName,      "read parameters from a file [default = 'scalingTest.xml']");
  bool        printTimings      = true;              clp.setOption("timings", "notimings",  &printTimings,     "print timings to screen");

  std::string mapFile;                               clp.setOption("map",                   &mapFile,          "map data file");
  std::string matrixFile;                            clp.setOption("matrix",                &matrixFile,       "matrix data file");
  int  optNraps   = 5;  clp.setOption("nraps",                &optNraps,   "number of RAPS to perform");
  bool optTimings = true; clp.setOption("timings", "notimings", &optTimings, "print timings to screen");
  bool optImplicitTranspose = true; clp.setOption("implicit", "explicit", &optImplicitTranspose, "whether to form R implicitly");
  bool optExport = false; clp.setOption("export", "noexport", &optExport, "write matrices to file");
  std::string optAggType="coupled"; clp.setOption("aggregation", &optAggType, "type of aggregation: \"uncoupled\",\"coupled\"");
  SC   optThreshold = 0.; clp.setOption("threshold", &optThreshold, "aggregation threshold");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  ParameterList paramList;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  bool isDriver = paramList.isSublist("Run1");
  if (isDriver) {
    // update galeriParameters with the values from the XML file
    ParameterList& realParams = galeriParameters.GetParameterList();

    for (ParameterList::ConstIterator it = realParams.begin(); it != realParams.end(); it++) {
      const std::string& name = realParams.name(it);
      if (paramList.isParameter(name))
        realParams.setEntry(name, paramList.getEntry(name));
    }
  }

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;
  comm->barrier();
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  RCP<Matrix>      A;
  RCP<const Map>   map;
  RCP<MultiVector> nullspace;
  RCP<MultiVector> coordinates;
  if (matrixFile.empty()) {
    galeriStream << "========================================================\n" << xpetraParameters << galeriParameters;

    // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a perfect square.
    // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.
    // If you don't want Galeri to do this, specify mx or my on the galeriList.
    std::string matrixType = galeriParameters.GetMatrixType();

    // Create map
    // In the future, we hope to be able to first create a Galeri problem, and then request map from it
    // At the moment, however, things are fragile as we hope that the Problem uses same map inside
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<int>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<int>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<int>("mz") << std::endl
                 << "========================================================" << std::endl;

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      galeriList.set("right boundary" , "Neumann");
      galeriList.set("bottom boundary", "Neumann");
      galeriList.set("top boundary"   , "Neumann");
      galeriList.set("front boundary" , "Neumann");
      galeriList.set("back boundary"  , "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" ||
        matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
    }

  } else {

    if (!mapFile.empty())
      map = Utils2::ReadMap(mapFile, xpetraParameters.GetLib(), comm);
    comm->barrier();

    if (lib == Xpetra::UseEpetra) {
      A = Utils::Read(matrixFile, map);

    } else {
      // Tpetra matrix reader is still broken, so instead we read in
      // a matrix in a binary format and then redistribute it
      const bool binaryFormat = true;
      A = Utils::Read(matrixFile, lib, comm, binaryFormat);

      RCP<Matrix> newMatrix = MatrixFactory::Build(map, 1);
      RCP<Import> importer  = ImportFactory::Build(A->getRowMap(), map);
      newMatrix->doImport(*A, *importer, Xpetra::INSERT);
      newMatrix->fillComplete();

      A.swap(newMatrix);
    }

    comm->barrier();

  }

  comm->barrier();
  tm = Teuchos::null;

  galeriStream << "Galeri complete.\n========================================================" << std::endl;

  const bool mustAlreadyExist = true;

  ParameterList mueluList, runList;
  if (isDriver) {
    runList   = paramList.sublist("Run1",  mustAlreadyExist);
    mueluList = runList  .sublist("MueLu", mustAlreadyExist);
  } else {
    mueluList = paramList;
  }

  if (nullspace.is_null()) {
    int blkSize = 1;
    if (mueluList.isSublist("Matrix")) {
      // Factory style parameter list
      const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
      if (operatorList.isParameter("PDE equations"))
        blkSize = operatorList.get<int>("PDE equations");

    } else if (paramList.isParameter("number of equations")) {
      // Easy style parameter list
      blkSize = paramList.get<int>("number of equations");
    }

    nullspace = MultiVectorFactory::Build(map, blkSize);
    for (int i = 0; i < blkSize; i++) {
      RCP<const Map> domainMap = A->getDomainMap();
      GO             indexBase = domainMap->getIndexBase();

      ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
      for (int j = 0; j < nsData.size(); j++) {
        GO GID = domainMap->getGlobalElement(j) - indexBase;

        if ((GID-i) % blkSize == 0)
          nsData[j] = Teuchos::ScalarTraits<SC>::one();
      }
    }
  }

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  out << galeriStream.str();

  Level fineLevel, coarseLevel;
  RAPFactory AcFact;
  AcFact.DisableMultipleCallCheck();
  RCP<SingleLevelFactoryBase> aggfact;
  RCP<TentativePFactory> ptentfact;
  RCP<CoalesceDropFactory> cdfact;
  RCP<SaPFactory>    PFact;
  RCP<TransPFactory> RFact;
  RCP<CoarseMapFactory> mapfact;

  {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("RAPScalingTest: 2 - Setup")));

    MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel); // set a default FactoryManager
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", coordinates);

    cdfact = rcp(new CoalesceDropFactory());
    ParameterList cdlist = *(cdfact->GetValidParameterList());
    cdlist.set("aggregation: drop tol",optThreshold);
    cdlist.set("aggregation: drop scheme","distance laplacian");
    cdlist.set("lightweight wrap",true);
    if (optAggType=="uncoupled") {
      aggfact = rcp(new UncoupledAggregationFactory());
    } else {
      aggfact = rcp(new CoupledAggregationFactory());
    }
    aggfact->SetFactory("Graph",cdfact);
    mapfact = rcp(new CoarseMapFactory() );
    mapfact->SetFactory("Aggregates",aggfact);
    ptentfact = rcp(new TentativePFactory());
    ptentfact->SetFactory("Aggregates",aggfact);
    ptentfact->SetFactory("CoarseMap",mapfact);
    cdfact->SetParameterList(cdlist);
    PFact = rcp(new SaPFactory());
    PFact->SetFactory("P",ptentfact);
    coarseLevel.Request("P", PFact.get());
    PFact->Build(fineLevel, coarseLevel);

    RFact = rcp(new TransPFactory());

    ParameterList Aclist = *(AcFact.GetValidParameterList());
    Aclist.set("transpose: use implicit", optImplicitTranspose);
    AcFact.SetParameterList(Aclist);

    AcFact.SetFactory("P", PFact);
    if (optImplicitTranspose==false)
      AcFact.SetFactory("R", RFact);
    tm = Teuchos::null;
  }

  RCP<Time> RAPKernelTimer = TimeMonitor::getNewTimer("RAPScalingTest: 3 - RAP kernel"); // re-use the same timer in the loop

  for (int i=0; i<optNraps; ++i) {
    coarseLevel.Request("A", &AcFact);
    {
      TimeMonitor tm2(*RAPKernelTimer);
      if (optImplicitTranspose==false) {
        RFact->SetFactory("P", PFact);
        coarseLevel.Request("R", RFact.get());
        RFact->Build(fineLevel, coarseLevel);
      }
      RCP<Matrix> Ac = coarseLevel.Get< RCP<Matrix> >("A", &AcFact);
    }
    coarseLevel.Release("A", &AcFact);
    if (optImplicitTranspose==false)
      coarseLevel.Release("R", RFact.get());
  }

  if (optExport) {
    std::string fileName = "A.m";
    Utils::Write(fileName,*A);
  }

  if (optTimings) {
    Teuchos::TableFormat &format = TimeMonitor::format();
    format.setPrecision(25);
    TimeMonitor::summarize();
  }

} //main
