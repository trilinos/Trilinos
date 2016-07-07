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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_VisualizationHelpers.hpp>
#include <MueLu_FacadeClassFactory.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ExportVTK : public MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    ExportVTK() {};

  public:

    void writeFile(std::ofstream& fout, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coordinates, Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& sol)
    {
      using namespace std;
      typedef MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> VH;
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > nodeMap = coordinates->getMap();
      std::vector<LocalOrdinal> vertices;
      std::vector<LocalOrdinal> geomSize;
      LocalOrdinal numFineNodes = Teuchos::as<LocalOrdinal>(coordinates->getLocalLength());

      vertices.reserve(numFineNodes);
      geomSize.reserve(numFineNodes);
      for(LocalOrdinal i = 0; i < numFineNodes; i++)
      {
        vertices.push_back(i);
        geomSize.push_back(1);
      }

      Teuchos::ArrayRCP<const double> xCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(0));
      Teuchos::ArrayRCP<const double> yCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(1));
      Teuchos::ArrayRCP<const double> zCoords = Teuchos::null;
      if(coordinates->getNumVectors() == 3) {
        zCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(2));
      }
      Teuchos::ArrayRCP<const double> solData = Teuchos::arcp_reinterpret_cast<const double>(sol->getData(0));

      std::vector<int> uniqueFine = this->makeUnique(vertices);
      this->writeFileVTKOpening(fout, uniqueFine, geomSize);
      this->writeFileVTKNodes(fout, uniqueFine, nodeMap);

      std::string indent = "      ";
      fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << std::endl;
      fout << indent;
      int myRank = coordinates->getMap()->getComm()->getRank();
      for(int i = 0; i < int(uniqueFine.size()); i++)
      {
        fout << myRank << " ";
        if(i % 20 == 19)
          fout << std::endl << indent;
      }
      fout << std::endl;
      fout << "        </DataArray>" << std::endl;
      // solution vector
      fout << "        <DataArray type=\"Float64\" Name=\"Solution\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      fout << indent;
      for(int i = 0; i < int(uniqueFine.size()); i++)
      {
        size_t numVec = coordinates->getNumVectors();
        for(int d=0; d<numVec; d++)
          fout << solData[i*numVec+d] << " ";
        if(i % 3 == 0)
          fout << std::endl << indent;
      }
      fout << std::endl;
      fout << "        </DataArray>" << std::endl;
      fout << "      </PointData>" << std::endl; // that is annoying

      this->writeFileVTKCoordinates(fout, uniqueFine, xCoords, yCoords, zCoords, coordinates->getNumVectors());
      this->writeFileVTKCells(fout, uniqueFine, vertices, geomSize);
      this->writeFileVTKClosing(fout);
      fout.close();

    };

  };

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& fancyout = *fancy;
    fancyout.setOutputToRootOnly(0);

    MueLu::FacadeClassFactory<SC,LO,GO,NO> test();

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

    //GO nx = 100, ny = 100, nz = 100;
    //Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    std::string xmlFileName = "driver.xml";      clp.setOption("xml",                   &xmlFileName,     "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");
    std::string belosFileName = "";              clp.setOption("belosXml",              &belosFileName,   "read parameters for Belos from a file.");
    int    amgAsPrecond     = 1;                 clp.setOption("precond",               &amgAsPrecond,     "apply multigrid as preconditioner");
    int    amgAsSolver      = 0;                 clp.setOption("fixPoint",              &amgAsSolver,      "apply multigrid as solver");
    bool   printTimings     = true;              clp.setOption("timings", "notimings",  &printTimings,     "print timings to screen");
    int    blockedSystem    = 0;                 clp.setOption("split",                 &blockedSystem,    "split system matrix into 2x2 system (default=0)");
    int    useThyraGIDs     = 0;                 clp.setOption("thyra",                 &useThyraGIDs,     "use Thyra style numbering of GIDs.");
    int    writeMatricesOPT = -2;                clp.setOption("write",                 &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");
    double tol              = 1e-6;             clp.setOption("tol",                   &tol,              "solver convergence tolerance");
    std::string krylovMethod = "gmres"; clp.setOption("krylov",                   &krylovMethod,     "outer Krylov method");
    int maxIts = 100; clp.setOption("maxits",           &maxIts,   "maximum number of Krylov iterations");
    int output = 1; clp.setOption("output",           &output,   "how often to print Krylov residual history");
    std::string matrixFileName = "crada1/crada_A.mm"; clp.setOption("matrixfile",            &matrixFileName,   "matrix market file containing matrix");
    std::string rhsFileName = "crada1/crada_b.mm";    clp.setOption("rhsfile",               &rhsFileName,      "matrix market file containing right-hand side");
    std::string nspFileName = "crada1/crada_ns.mm";   clp.setOption("nspfile",               &nspFileName,      "matrix market file containing fine level null space");
    std::string cooFileName = "crada1/crada_coordinates.mm"; clp.setOption("coordinatesfile",&cooFileName,      "matrix market file containing fine level coordinates");
    std::string spcFileName = "crada1/crada_special.mm"; clp.setOption("specialfile",        &spcFileName,      "matrix market file containing fine level special dofs");
    int nPDE = 3; clp.setOption("numpdes",           &nPDE,   "number of PDE equations");
    int nNspVectors = 3; clp.setOption("numnsp", &nNspVectors, "number of nullspace vectors. Only used if null space is read from file. Must be smaller or equal than the number of null space vectors read in from file.");
    std::string convType = "r0"; clp.setOption("convtype",       &convType,         "convergence type (r0 or none)");
    std::string strOutputFilename = ""; clp.setOption("output",  &strOutputFilename,"filename prefix for output file name. If empty, no output is written.");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;

    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));

    RCP<Matrix> A = Teuchos::null;
    if (matrixFileName != "") {
      fancyout << "Read matrix from file " << matrixFileName << std::endl;
      RCP<Matrix> Atest = Xpetra::IO<SC,LO,GO,Node>::Read(std::string(matrixFileName), xpetraParameters.GetLib(), comm);
      RCP<const Map>   maptest = Atest->getRowMap();

      // re-read matrix and make sure it is properly distributed over all processors
      // make sure that all DOF rows per node are located on same processor.
      GO numGlobalNodes = maptest->getGlobalNumElements() / nPDE;
      GO numMyNodes     = numGlobalNodes / comm->getSize();
      if (comm->getRank() == comm->getSize()-1)
        numMyNodes = numGlobalNodes - (comm->getSize()-1) * numMyNodes;
      LO numMyDofs = Teuchos::as<LocalOrdinal>(numMyNodes) * nPDE;

      // construct new row map for equally distributing the matrix rows without
      // splitting dofs that belong to the same node
      Teuchos::Array<GO> myDofGIDs(numMyDofs);
      for(LO i = 0; i < numMyDofs; i++) {
        if (comm->getRank() == comm->getSize()-1)
          myDofGIDs[i] = Teuchos::as<GO>(maptest->getGlobalNumElements()) - Teuchos::as<GO>(numMyDofs - i);
        else
          myDofGIDs[i] = Teuchos::as<GO>(comm->getRank() * numMyDofs) + Teuchos::as<GO>(i);
      }
      RCP<const Map> Arowmap = MapFactory::Build (xpetraParameters.GetLib(),Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),myDofGIDs(),0,comm);

      A = Xpetra::IO<SC,LO,GO,Node>::Read(matrixFileName,Arowmap);
      A->SetFixedBlockSize(Teuchos::as<LocalOrdinal>(nPDE));
    }
    RCP<const Map>   map = A->getRowMap();
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(),nPDE);
    A->SetFixedBlockSize(nPDE);
    fancyout << "#pdes = " << A->GetFixedBlockSize() << std::endl;

    if (nspFileName != "") {
      fancyout << "Read null space from file " << nspFileName << std::endl;
      nullspace = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(std::string(nspFileName), A->getRowMap());
      fancyout << "Found " << nullspace->getNumVectors() << " null space vectors" << std::endl;
      if (nNspVectors > Teuchos::as<int>(nullspace->getNumVectors())) {
        fancyout << "Set number of null space vectors from " << nNspVectors << " to " << nullspace->getNumVectors() << " as only " << nullspace->getNumVectors() << " are provided by " << nspFileName << std::endl;
        nNspVectors = nullspace->getNumVectors();
      }
      if (nNspVectors < 1) {
        fancyout << "Set number of null space vectors from " << nNspVectors << " to " << nullspace->getNumVectors() << ". Note: we need at least one null space vector!!!" << std::endl;
        nNspVectors = nullspace->getNumVectors();
      }
      if (nNspVectors < Teuchos::as<int>(nullspace->getNumVectors())) {
        RCP<MultiVector> temp = MultiVectorFactory::Build(A->getDomainMap(),nNspVectors);
        for(int j=0; j<nNspVectors; j++) {
          Teuchos::ArrayRCP<SC> tempData = temp->getDataNonConst(j);
          Teuchos::ArrayRCP<const SC> nsData   = nullspace->getData(j);
          for (int i=0; i<nsData.size(); ++i) {
            tempData[i] = nsData[i];
          }
        }
        nullspace = Teuchos::null;
        nullspace = temp;
      }
    } else {
      if (nPDE == 1)
        nullspace->putScalar( Teuchos::ScalarTraits<SC>::one() );
      else {
        for (int i=0; i<nPDE; ++i) {
          Teuchos::ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
          for (int j=0; j<nsData.size(); ++j) {
            GO gel = A->getDomainMap()->getGlobalElement(j) - A->getDomainMap()->getIndexBase();
            if ((gel-i) % nPDE == 0)
              nsData[j] = Teuchos::ScalarTraits<SC>::one();
          }
        }
      }
    }

    RCP<MultiVector> coordinates = Teuchos::null; //MultiVectorFactory::Build(A->getDomainMap(),1);
    if (cooFileName != "") {
      std::vector<GO> myGIDs (map->getNodeNumElements() / A->GetFixedBlockSize());
      // reconstruct map for coordinates
      for(LO r = 0; r < Teuchos::as<LO>(map->getNodeNumElements() / A->GetFixedBlockSize()); ++r) {
        GO gid = map->getGlobalElement(r * A->GetFixedBlockSize());
        myGIDs[r] = gid;
      }

      GO gCntGIDs  = 0;
      GO glCntGIDs = Teuchos::as<GlobalOrdinal>(myGIDs.size());
      MueLu_sumAll(comm,glCntGIDs,gCntGIDs);

      Teuchos::Array<GlobalOrdinal> eltList(myGIDs);
      RCP<const Map> myCoordMap = MapFactory::Build (xpetraParameters.GetLib(),gCntGIDs,eltList(),0,comm);

      fancyout << "Read fine level coordinates from file " << cooFileName << std::endl;
      coordinates = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(std::string(cooFileName), myCoordMap);
      fancyout << "Found " << nullspace->getNumVectors() << " coordinate vectors of length " << myCoordMap->getGlobalNumElements() << std::endl;
    }

    // shouldn't these be const?
    RCP<Map> mySpecialMap = Teuchos::null;
    if (spcFileName != "") {
      // read file on each processor and pick out the special dof numbers which belong to the current proc
      std::ifstream infile(spcFileName);
      std::string line;
      Teuchos::Array<GlobalOrdinal> mySpecialGids;
      Teuchos::Array<GlobalOrdinal> nonSpecialGids;
      GlobalOrdinal cnt = 0;   // count overall number of gids
      GlobalOrdinal mycnt = 0; // count only local gids
      while ( std::getline(infile, line)) {
        if(0 == line.find("%")) continue;
        if(0 == line.find(" ")) {
          cnt++;
          GlobalOrdinal gid;
          std::istringstream iss(line);
          iss >> gid;
          gid--; // note, that the matlab vector starts counting at 1 and not 0!
          if(map->isNodeGlobalElement(gid)) {
            mySpecialGids.push_back(gid);
            mycnt++;
          }
        }
      }

      //Teuchos::Array<GlobalOrdinal> eltList(mySpecialGids);
      mySpecialMap    = MapFactory::Build (xpetraParameters.GetLib(),cnt,mySpecialGids(),0,comm);

      // empty processors
      std::vector<size_t> lelePerProc(comm->getSize(),0);
      std::vector<size_t> gelePerProc(comm->getSize(),0);
      lelePerProc[comm->getRank()] = mySpecialMap->getNodeNumElements();
      Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&lelePerProc[0],&gelePerProc[0]);
      if(comm->getRank() == 0) {
        fancyout << "Distribution of " << cnt << " special dofs over processors:" << std::endl;
        fancyout << "Proc   #DOFs" << std::endl;
        for(int i=0; i<comm->getSize(); i++) {
         fancyout  << i << "      " << gelePerProc[i] << std::endl;
        }
      }
    }

    comm->barrier();
    tm = Teuchos::null;

    comm->barrier();
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1.2 - LHS and RHS initialization")));

    RCP<Vector> X = VectorFactory::Build(map,1);
    RCP<MultiVector> B = VectorFactory::Build(map,1);

    if (rhsFileName != "")
      B = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(std::string(rhsFileName), A->getRowMap());
    else
    {
      // we set seed for reproducibility
      X->setSeed(846930886);
      bool useSameRandomGen = false;
      X->randomize(useSameRandomGen);
      A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
      B->norm2(norms);
      //B->scale(1.0/norms[0]);
    }

    X->putScalar(zero);

    comm->barrier();
    tm = Teuchos::null;

    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1.5 - MueLu read XML and feed MueLu")));
    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);


    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

    // By default, we use Extreme. However, typically the xml file contains verbosity parameter
    // which is used instead
    H->SetDefaultVerbLevel(MueLu::Extreme);

    // =========================================================================
    // Prepare linear system
    // =========================================================================
    // We have the following different formulations:
    // - split the linear system in a 2x2 multiphysics problem using Xpetra style gids
    // - split the linear system in a 2x2 multiphysics problem using Thyra style gids
    // - solve the problem as a monolithik linear system
    if(blockedSystem == 1) {
      // split matrix and vectors

      // create map extractor
      Teuchos::Array<GlobalOrdinal> nonSpecialGids;
      for (size_t i = 0; i < map->getNodeNumElements(); i++) {
        GlobalOrdinal gid = map->getGlobalElement(i);
        if (mySpecialMap->isNodeGlobalElement(gid) == false) {
          nonSpecialGids.push_back(gid);
        }
      }
      //MapFactory::Build (xpetraParameters.GetLib(),Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),nonSpecialGids(),0,comm);
      std::vector<size_t> strInfo(1,nPDE);
      RCP<const Map> myStridedNonSpecialMap = StridedMapFactory::Build(xpetraParameters.GetLib(),Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),nonSpecialGids(),0,strInfo,comm);
      RCP<const Map> myStridedSpecialMap    = StridedMapFactory::Build(mySpecialMap, strInfo);
      //std::cout << "map " << map->getMaxAllGlobalIndex() << "nonspecial " << myNonSpecialMap->getMaxAllGlobalIndex() << " " << mySpecialMap->getMaxAllGlobalIndex() << std::endl;
      //std::cout << Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<int, Node> >(myNonSpecialMap)->getEpetra_Map() << std::endl;
      //std::cout << Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<int, Node> >(mySpecialMap)->getEpetra_Map() << std::endl;

      std::vector<Teuchos::RCP<const Map> > xmaps;
      xmaps.push_back(myStridedNonSpecialMap);
      xmaps.push_back(myStridedSpecialMap);

      Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map,xmaps);

      // split null space vectors
      RCP<MultiVector> nullspace1 = map_extractor->ExtractVector(nullspace,0);
      RCP<MultiVector> nullspace2 = map_extractor->ExtractVector(nullspace,1);

      bool bThyraMode = (useThyraGIDs==1) ? true : false;
      Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bOp =
          Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplitMatrix(*A,map_extractor,map_extractor,Teuchos::null,bThyraMode);

      // TODO plausibility checks
      // TODO set number of dofs per node

      // use bOp as A
      A = bOp;

      if(bThyraMode == false) {
        // use Xpetra style GIDs
        H->GetLevel(0)->Set("A",            Teuchos::rcp_dynamic_cast<Matrix>(bOp));
        H->GetLevel(0)->Set("Nullspace1",   nullspace1);
        H->GetLevel(0)->Set("Nullspace2",   nullspace2);
        H->GetLevel(0)->Set("Coordinates", coordinates); // TODO split coordinates for rebalancing! (or provide the full vector in the right map and split it in the factories!)
        if(mySpecialMap!=Teuchos::null) H->GetLevel(0)->Set("map SpecialMap", mySpecialMap);
      } else {
        // use Thyra style GIDs
        RCP<MultiVector> nsp1shrinked   = Xpetra::MatrixUtils<Scalar,LocalOrdinal, GlobalOrdinal, Node>::xpetraGidNumbering2ThyraGidNumbering(*nullspace1);
        RCP<MultiVector> nsp2shrinked   = Xpetra::MatrixUtils<Scalar,LocalOrdinal, GlobalOrdinal, Node>::xpetraGidNumbering2ThyraGidNumbering(*nullspace2);
        RCP<MultiVector> coordsshrinked = Xpetra::MatrixUtils<Scalar,LocalOrdinal, GlobalOrdinal, Node>::xpetraGidNumbering2ThyraGidNumbering(*coordinates);
        H->GetLevel(0)->Set("A",           Teuchos::rcp_dynamic_cast<Matrix>(bOp));
        H->GetLevel(0)->Set("Nullspace1",  nsp1shrinked);
        H->GetLevel(0)->Set("Nullspace2",  nsp2shrinked);
        H->GetLevel(0)->Set("Coordinates", coordsshrinked);  // TODO split coordinates for rebalancing! (or provide the full vector in the right map and split it in the factories!)
        if(mySpecialMap!=Teuchos::null) {
          RCP<const Map> specialmapshrinked = Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node>::shrinkMapGIDs(*mySpecialMap,*mySpecialMap);
          H->GetLevel(0)->Set("map SpecialMap", Teuchos::rcp_const_cast<Map>(specialmapshrinked));
        }

        // rearrange contents of rhs vector B
        RCP<const MapExtractor> thy_map_extractor = bOp->getRangeMapExtractor();
        RCP<MultiVector> Bshrinked = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(thy_map_extractor->getFullMap(),B->getNumVectors(),true);

        size_t numMaps = map_extractor->NumMaps();
        for(size_t k = 0; k < numMaps; k++) {
          RCP<MultiVector> partVec         = map_extractor->ExtractVector(B, k);
          RCP<MultiVector> newPartVec      = thy_map_extractor->getVector(k, B->getNumVectors(), false); // we need real GIDs (not zero based Thyra GIDs)
          // copy data
          for (size_t c = 0; c < partVec->getNumVectors(); c++) {
            Teuchos::ArrayRCP< const Scalar > data = partVec->getData(c);
            for (size_t r = 0; r < partVec->getLocalLength(); r++) {
              newPartVec->replaceLocalValue(Teuchos::as<LocalOrdinal>(r), c, data[r]);
            }
          }
          thy_map_extractor->InsertVector(*newPartVec, k, *Bshrinked, false);
        }

        B = Bshrinked;
        X = VectorFactory::Build(thy_map_extractor->getFullMap(),1);
        X->putScalar(zero);

        // TODO the ordering of the solution vector X is different!
      }
    } else {
      // standard (non-blocked) case
      H->GetLevel(0)->Set("A",           A);
      H->GetLevel(0)->Set("Nullspace",   nullspace);
      H->GetLevel(0)->Set("Coordinates", coordinates);
      if(mySpecialMap!=Teuchos::null) H->GetLevel(0)->Set("map SpecialMap", mySpecialMap);
    }

    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - MueLu Setup")));

    mueLuFactory.SetupHierarchy(*H);

    comm->barrier();
    tm = Teuchos::null;

    // Print out the hierarchy stats. We should not need this line, but for some reason the
    // print out in the hierarchy construction does not work.
    H->print(fancyout);

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================

    tm = Teuchos::null;

    if (writeMatricesOPT > -2)
      H->Write(writeMatricesOPT, writeMatricesOPT);

    comm->barrier();
    if (amgAsSolver) {
      tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Fixed Point Solve")));

      H->IsPreconditioner(false);
      Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
      norms = Utilities::ResidualNorm(*A,*X,*B);
      std::cout << "                iter:    0           residual = " << norms[0] << std::endl;
      for (int i=0; i< maxIts; ++i) {
        H->Iterate(*B, *X);
        norms = Utilities::ResidualNorm(*A,*X,*B);
        std::cout << "                iter:    " << i+1 << "           residual = " << norms[0] << std::endl;
      }

    } else if (amgAsPrecond) {
#ifdef HAVE_MUELU_BELOS
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Belos Solve")));
      // Operator and Multivector type that will be used with Belos
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;
      H->IsPreconditioner(true);

      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
      Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator
      // Construct a Belos LinearProblem object
      RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
      belosProblem->setRightPrec(belosPrec);
      bool set = belosProblem->setProblem();
      if (set == false) {
        fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos parameter list
      Teuchos::ParameterList belosList;
      if(belosFileName == "") {
        belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
        belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
        belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList.set("Output Frequency",      output);
        belosList.set("Output Style",          Belos::Brief);
        if (convType == "none") {
          belosList.set("Explicit Residual Scaling",  "None");
          belosList.set("Implicit Residual Scaling",  "None");
        }
      } else {
        Teuchos::updateParametersFromXmlFileAndBroadcast(belosFileName, Teuchos::Ptr<Teuchos::ParameterList>(&belosList), *comm);

        /*belosList.set("Output Style",          Belos::User);
        belosList.set("User Status Tests Combo Type", "OR");
        Teuchos::ParameterList& userList = belosList.sublist("User Status Tests");
        userList.set("Test Type","Combo");
        userList.set("Combo Type","AND");
        userList.set("Number of Tests",3);
        Teuchos::ParameterList& userList1 = userList.sublist("Test 0");
        userList1.set("Tag","1 FULL");
        userList1.set("Test Type","ResidualNorm");
        userList1.set("Convergence Tolerance",tol);
        userList1.set("Scaling Type","None");
        userList1.set("Residual Norm","OneNorm");
        userList1.set("Scaling Norm","OneNorm");
        Teuchos::ParameterList& userList2 = userList.sublist("Test 1");
        userList2.set("Tag","2 PRIM");
        userList2.set("Test Type","PartialResidualNorm");
        userList2.set("Block index",0);
        userList2.set("Convergence Tolerance",tol);
        userList2.set("Scaling Type","None");
        userList2.set("Residual Norm","OneNorm");
        userList2.set("Scaling Norm","OneNorm");
        Teuchos::ParameterList& userList3 = userList.sublist("Test 2");
        userList3.set("Tag","3 SECOND");
        userList3.set("Test Type","PartialResidualNorm");
        userList3.set("Block index",1);
        userList3.set("Convergence Tolerance",tol);
        userList3.set("Scaling Type","None");
        userList3.set("Residual Norm","OneNorm");
        userList3.set("Scaling Norm","OneNorm");*/
      }

      // Create an iterative solver manager
      RCP< Belos::SolverManager<SC, MV, OP> > solver;
      if (krylovMethod == "cg") {
        solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
      } else if (krylovMethod == "gmres") {
        solver = rcp(new Belos::PseudoBlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid Krylov method.  Options are \"cg\" or \" gmres\".");
      }

      // Perform solve
      Belos::ReturnType ret = Belos::Unconverged;
      try {
        ret = solver->solve();

        // Get the number of iterations for this solve.
        fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      } catch(...) {
        fancyout << std::endl << "ERROR:  Belos threw an error! " << std::endl;
      }
      // Check convergence
      if (ret != Belos::Converged)
        fancyout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
      else
        fancyout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
#endif //ifdef HAVE_MUELU_BELOS
    }
    comm->barrier();
    tm = Teuchos::null;

    if (strOutputFilename.empty() == false) {
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Export solution in VTK")));
      if(comm->getSize() > 1) {
        std::stringstream ss;
        ss << "-proc" << comm->getRank();
        strOutputFilename.append(ss.str());
      }
      strOutputFilename.append(".vtu");

      std::ofstream fout(strOutputFilename);
      ExportVTK<Scalar,LocalOrdinal,GlobalOrdinal,Node> expVTK;
      expVTK.writeFile(fout,coordinates,X);
      fout.close();
    }

    // print timings

    comm->barrier();
    tm = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
      MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main
