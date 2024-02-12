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
#include <cstdio>
#include <unistd.h>
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

bool compare_to_gold_all_ranks(int myRank, const std::string& baseFile) {
  bool failed = false;

  // Create a copy of outputs
  std::string cmd = "cp -f ";
  int retVal;
  retVal = system((cmd + baseFile + ".gold " + baseFile + ".gold_filtered").c_str());
  TEUCHOS_TEST_FOR_EXCEPTION(retVal != 0, std::runtime_error, "goldfile copy failed!");
  // system((cmd + baseFile + ".out "  + baseFile + ".out_filtered").c_str());
  retVal = system((cmd + baseFile + ".vtu " + baseFile + ".out_filtered").c_str());
  TEUCHOS_TEST_FOR_EXCEPTION(retVal != 0, std::runtime_error, "goldfile copy failed!");

  // Run comparison (ignoring whitespaces)
  cmd     = "diff -u -w -I\"^\\s*$\" " + baseFile + ".gold_filtered " + baseFile + ".out_filtered";
  int ret = system(cmd.c_str());
  std::cout << ret << std::endl;
  if (ret)
    failed = true;

  std::cout << baseFile << ": " << (ret ? "failed" : "passed") << std::endl;

  return !failed;
}

std::string replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) {
  while (1) {
    const int pos = result.find(replaceWhat);
    if (pos == -1)
      break;
    result.replace(pos, replaceWhat.size(), replaceWithWhat);
  }
  return result;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  bool success  = false;
  bool aggMatch = false;
  bool verbose  = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    typedef Teuchos::ScalarTraits<SC> STS;
    SC one = STS::one();

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

    std::string xmlFileName = "viztest.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file");

    std::string mapFile;
    clp.setOption("map", &mapFile, "map data file");
    std::string colMapFile;
    clp.setOption("colmap", &colMapFile, "colmap data file");
    std::string domainMapFile;
    clp.setOption("domainmap", &domainMapFile, "domainmap data file");
    std::string rangeMapFile;
    clp.setOption("rangemap", &rangeMapFile, "rangemap data file");
    std::string matrixFile;
    clp.setOption("matrix", &matrixFile, "matrix data file");
    std::string coordFile;
    clp.setOption("coords", &coordFile, "coordinates data file");
    std::string nullFile;
    clp.setOption("nullspace", &nullFile, "nullspace data file");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);

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

    RCP<Matrix> A;
    RCP<const Map> map;
    RCP<MultiVector> coordinates;
    RCP<MultiVector> nullspace;
    int ndims = 2;
    if (matrixFile.empty()) {
      galeriStream << "========================================================\n"
                   << xpetraParameters << galeriParameters;

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

      // Create map and coordinates
      // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
      // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
      if (matrixType == "Laplace1D") {
        map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, galeriList);

      } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
                 matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
        map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", map, galeriList);

      } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
        map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("3D", map, galeriList);
        ndims       = 3;
      }

      // Expand map to do multiple DOF per node for block problems
      if (matrixType == "Elasticity2D")
        map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
      if (matrixType == "Elasticity3D")
        map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

      galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
                   << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
                   << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
                   << "========================================================" << std::endl;

      if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
        // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
        galeriList.set("right boundary", "Neumann");
        galeriList.set("bottom boundary", "Neumann");
        galeriList.set("top boundary", "Neumann");
        galeriList.set("front boundary", "Neumann");
        galeriList.set("back boundary", "Neumann");
      }

      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
      A = Pr->BuildMatrix();

      if (matrixType == "Elasticity2D" ||
          matrixType == "Elasticity3D") {
        nullspace = Pr->BuildNullspace();
        A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
      }

    } else {
      if (!mapFile.empty())
        map = Xpetra::IO<SC, LO, GO, Node>::ReadMap(mapFile, lib, comm);
      comm->barrier();

      const bool binaryFormat = false;

      if (!binaryFormat && !map.is_null()) {
        RCP<const Map> colMap    = (!colMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(colMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> domainMap = (!domainMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(domainMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> rangeMap  = (!rangeMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(rangeMapFile, lib, comm) : Teuchos::null);
        A                        = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, map, colMap, domainMap, rangeMap);

      } else {
        A = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, lib, comm, binaryFormat);

        if (!map.is_null()) {
          RCP<Matrix> newMatrix = MatrixFactory::Build(map, 1);
          RCP<Import> importer  = ImportFactory::Build(A->getRowMap(), map);
          newMatrix->doImport(*A, *importer, Xpetra::INSERT);
          newMatrix->fillComplete();

          A.swap(newMatrix);
        }
      }
      map = A->getMap();

      comm->barrier();

      if (!coordFile.empty()) {
        // NOTE: currently we only allow reading scalar matrices, thus coordinate
        // map is same as matrix map
        coordinates = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(coordFile, map);
      }

      if (!nullFile.empty())
        nullspace = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(nullFile, map);
    }

    comm->barrier();
    tm = Teuchos::null;

    galeriStream << "Galeri complete.\n========================================================" << std::endl;

    std::string matrixType = galeriParameters.GetMatrixType();
    std::string aggVizType = paramList.get<std::string>("aggregation: output file: agg style");
    aggVizType.erase(std::remove_if(aggVizType.begin(), aggVizType.end(), ::isspace), aggVizType.end());
    if (ndims == 2)
      paramList.set<std::string>("aggregation: output filename", "MPI-Viz-Output-2D-Level%LEVELID-Proc%PROCID");
    else if (ndims == 3) {
      if (comm->getSize() > 1)
        paramList.set<std::string>("aggregation: output filename", "MPI-Viz-Output-" + matrixType + "-" + aggVizType + "-Level%LEVELID-Proc%PROCID");
      else
        paramList.set<std::string>("aggregation: output filename", "MPI-Viz-Output-" + matrixType + "-" + aggVizType + "-Level%LEVELID");
    }

    if (nullspace.is_null()) {
      int blkSize = 1;
      if (paramList.isSublist("Matrix")) {
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
        GO indexBase             = domainMap->getIndexBase();

        ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
        for (int j = 0; j < nsData.size(); j++) {
          GO GID = domainMap->getGlobalElement(j) - indexBase;

          if ((GID - i) % blkSize == 0)
            nsData[j] = Teuchos::ScalarTraits<SC>::one();
        }
      }
    }

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out       = *fancy;
    out.setOutputToRootOnly(0);

    out << galeriStream.str();

    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    comm->barrier();
    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));

    RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(paramList));

    comm->barrier();
    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

    RCP<Hierarchy> H;
    A->SetMaxEigenvalueEstimate(-one);

    H = mueLuFactory->CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    if (!coordinates.is_null())
      H->GetLevel(0)->Set("Coordinates", coordinates);
    mueLuFactory->SetupHierarchy(*H);

    comm->barrier();
    tm = Teuchos::null;

    // =========================================================================
    // Check aggs
    // =========================================================================
    std::string filenameToWrite;
    if (comm->getSize() > 1)
      filenameToWrite = "MPI-Viz-Output-" + matrixType + "-" + aggVizType + "-Level0-Proc%PROCID";
    else
      filenameToWrite = "MPI-Viz-Output-" + matrixType + "-" + aggVizType + "-Level0";
    std::string outfileName = replaceAll(filenameToWrite, "%PROCID", MueLu::toString(comm->getRank()));
    aggMatch                = compare_to_gold_all_ranks(comm->getRank(), outfileName);

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    comm->barrier();

    success = aggMatch;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
