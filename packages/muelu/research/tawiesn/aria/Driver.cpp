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
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
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

#include <MueLu_MutuallyExclusiveTime.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

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
    Teuchos::FancyOStream& out = *fancy;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================

    //GO nx = 100, ny = 100, nz = 100;
    //Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    std::string xmlFileName = "driver.xml";      clp.setOption("xml",                   &xmlFileName,     "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    // =========================================================================
    // Read in problem
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;


    GlobalOrdinal nGlobalNodes = 0;
    if(comm->getRank() == 0) {
      std::ifstream data_file;
      data_file.open("xxx.mm");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file.good() == false, MueLu::Exceptions::RuntimeError,"Problem opening file xxx.mm");
      std::string line;
      getline(data_file,line);
      data_file >> nGlobalNodes;
      data_file.close();
    }
    Teuchos::broadcast(*comm,0,1,&nGlobalNodes);

    out << "Found " << nGlobalNodes << " nodes " << std::endl;

    // each processor reads in its dofPresent array
    LocalOrdinal nNodes = 0;
    LocalOrdinal nDofs  = 0;
    int maxDofPerNode = -1;
    Teuchos::ArrayRCP<LocalOrdinal> dofPresent;
    {
      FILE* data_file;
      std::stringstream ss; ss << comm->getSize() << "dofPresent" << comm->getRank();
      data_file = fopen(ss.str().c_str(),"r");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file == NULL, MueLu::Exceptions::RuntimeError,"Problem opening file " << ss.str());
      // read in first line containing number of local nodes and maxDofPerNode
      fscanf(data_file,"%d %d\n", &nNodes, &maxDofPerNode);
      dofPresent = Teuchos::ArrayRCP<LocalOrdinal>(nNodes * maxDofPerNode,0);
      // loop over all local nodes
      for(LocalOrdinal i = 0; i < nNodes; i++) {
        for(int j=0; j<maxDofPerNode; j++) {
          int tmp = -1;
          fscanf(data_file,"%d",&tmp);
          if(tmp == 1) {
            dofPresent[i*maxDofPerNode+j] = 1; nDofs++;
          } else {
            dofPresent[i*maxDofPerNode+j] = 0;
          }
        }
      }
      fclose(data_file);
    }

    out << "PROC " << comm->getRank() << "/" << comm->getSize() << " " << nNodes << " " << maxDofPerNode << std::endl;


    /*for(LocalOrdinal i = 0; i < nNodes; i++) {
      for(int j=0; j<maxDofPerNode; j++) {
        std::cout << i << " " << dofPresent[i*maxDofPerNode+j] << std::endl;
      }
    }*/


    // read in dofs

    Teuchos::Array<GlobalOrdinal> dofGlobals;
    Teuchos::Array<GlobalOrdinal> nodalGlobals;  // nodal GIDs for laplacian (with holes)
    {
      FILE* data_file;
      std::stringstream ss; ss << comm->getSize() << "proc" << comm->getRank();
      data_file = fopen(ss.str().c_str(),"r");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file == NULL, MueLu::Exceptions::RuntimeError,"Problem opening file " << ss.str());
      // read in first line containing number of local nodes and maxDofPerNode
      LocalOrdinal tmpDofs = 0;
      fscanf(data_file,"%d\n", &tmpDofs);
      TEUCHOS_TEST_FOR_EXCEPTION(tmpDofs != nDofs, MueLu::Exceptions::RuntimeError,"Number of gid entries in map file is " << tmpDofs << " but should be " << nDofs);

      dofGlobals = Teuchos::Array<GlobalOrdinal>(nDofs);

      // loop over all local nodes
      for(GlobalOrdinal i = 0; i < nDofs; i++) {
        int data;
        fscanf(data_file,"%d",&data);
        dofGlobals[i] = Teuchos::as<GlobalOrdinal>(data);
      }
      fclose(data_file);

      nodalGlobals = Teuchos::Array<GlobalOrdinal>(nNodes);

      GlobalOrdinal count = nDofs - 1;
      for(GlobalOrdinal i = nNodes - 1; i>=0; i--) {
        for(int j = maxDofPerNode-1; j >=0; j--) {
          if(dofPresent[i*maxDofPerNode+j] == 1)
            nodalGlobals[i] = dofGlobals[count--];
        }
      }
    }

    Teuchos::RCP<Map> LapMap = MapFactory::Build(lib,nGlobalNodes,nodalGlobals(),0,comm);

    //LapMap->describe(out, Teuchos::VERB_EXTREME);

    //

    {
      FILE* data_file;
      std::stringstream ss; ss << comm->getSize() << "ProcLinear" << comm->getRank();
      data_file = fopen(ss.str().c_str(),"r");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file == NULL, MueLu::Exceptions::RuntimeError,"Problem opening file " << ss.str());
      // loop over all local nodes
      for(GlobalOrdinal i = 0; i < nNodes; i++) {
        int data;
        fscanf(data_file,"%d",&data);
        nodalGlobals[i] = Teuchos::as<GlobalOrdinal>(data);
      }
      fclose(data_file);
      for(GlobalOrdinal i = 0; i < nNodes; i++) {
        nodalGlobals[i] = nodalGlobals[i] - 1;
      }
    }

    Teuchos::RCP<Map> LinearMap = MapFactory::Build(lib,nGlobalNodes,nodalGlobals(),0,comm);
    //LinearMap->describe(out, Teuchos::VERB_EXTREME);LapMap->describe(out, Teuchos::VERB_EXTREME);

    //LapMap->describe(out, Teuchos::VERB_EXTREME);LapMap->describe(out, Teuchos::VERB_EXTREME);

    Teuchos::RCP<Import> Importer = ImportFactory::Build(LinearMap, LapMap);

    Teuchos::RCP<MultiVector> xpetraXXX = MultiVectorFactory::Build(LapMap,1);
    Teuchos::RCP<MultiVector> xpetraYYY = MultiVectorFactory::Build(LapMap,1);

    RCP<MultiVector> temp = IO::ReadMultiVector ("xxx.mm", LinearMap);
    xpetraXXX->doImport(*temp, *Importer, Xpetra::INSERT);
    temp = IO::ReadMultiVector ("yyy.mm", LinearMap);
    xpetraYYY->doImport(*temp, *Importer, Xpetra::INSERT);

    Teuchos::RCP<MultiVector> coordinates = MultiVectorFactory::Build(LapMap,2);
    Teuchos::ArrayRCP< const Scalar > srcX = xpetraXXX->getData(0);
    Teuchos::ArrayRCP< const Scalar > srcY = xpetraYYY->getData(0);
    Teuchos::ArrayRCP< Scalar > dataX = coordinates->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > dataY = coordinates->getDataNonConst(1);
    for(decltype(coordinates->getLocalLength()) i = 0; i < coordinates->getLocalLength(); i++) {
      dataX[i] = srcX[i];
      dataY[i] = srcY[i];
    }

    // read in matrix
    GlobalOrdinal NumRows = 0; // number of rows in matrix
    GlobalOrdinal NumElements = 0; // number of elements in matrix
    GlobalOrdinal Offset = 0; // number of elements in matrix
    Teuchos::RCP<Matrix> SerialMatrix = Teuchos::null;
    {
      std::ifstream data_file;
      if (comm->getRank() == 0) {
        // proc 0 reads the number of rows, columns and nonzero elements
        data_file.open("theMatrix");
        TEUCHOS_TEST_FOR_EXCEPTION(data_file.good() == false, MueLu::Exceptions::RuntimeError,"Problem opening file theMatrix");
        data_file >> NumRows;
        data_file >> NumElements;
        data_file >> Offset;

      }

      // create a standard (serial) row map for the Matrix data (only elements on proc 0)
      Teuchos::RCP<Map> SerialMap = MapFactory::Build(lib,NumRows,NumRows,0,comm);
      SerialMatrix = MatrixFactory::Build(SerialMap, 33);
      if (comm->getRank() == 0) {
        for(decltype(NumElements) i = 0; i < NumElements; ++i) {
          GlobalOrdinal row;
          GlobalOrdinal col;
          Scalar val;
          data_file >> row;
          data_file >> col;
          data_file >> val;
          row -= Offset;
          col -= Offset;
          SerialMatrix->insertGlobalValues(row,Teuchos::Array<GlobalOrdinal>(1,col)(),Teuchos::Array<Scalar>(1,val)());
          if (i % 10000 == 0) {
            double percent = i * 100 / NumElements;
            out << "Read matrix: " << percent << "% complete" << std::endl;
          }
        }
        SerialMatrix->fillComplete();

        data_file.close();
      }

      //SerialMatrix->describe(out, Teuchos::VERB_EXTREME);

    } // end read in file

    // distribute map and matrix over processors
    Teuchos::RCP<const Map>    DistributedMap    = Teuchos::null;
    Teuchos::RCP<Matrix> DistributedMatrix = Teuchos::null;
    {
      if(comm->getSize() > 1) {
        Teuchos::broadcast(*comm,0,1,&NumRows);

        DistributedMap    = MapFactory::Build(lib,NumRows,dofGlobals(),0,comm);
        DistributedMatrix = MatrixFactory::Build(DistributedMap,33);

        Teuchos::RCP<Import> dofMatrixImporter = ImportFactory::Build(SerialMatrix->getRowMap(), DistributedMap);

        DistributedMatrix->doImport(*SerialMatrix,*dofMatrixImporter,Xpetra::INSERT);
        DistributedMatrix->fillComplete();
      } else {
        DistributedMap = SerialMatrix->getRowMap();
        DistributedMatrix = SerialMatrix;
      }

    } // end distribute matrix

    // read in global vectors (e.g. rhs)
    GlobalOrdinal nGlobalDof = 0;
    GlobalOrdinal nLocalDofs = Teuchos::as<GlobalOrdinal>(nDofs);

    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,comm->getSize(),&nLocalDofs,&nGlobalDof);

    Teuchos::RCP<const Map> dofLinearMap = Teuchos::null;
    {
      FILE* data_file;
      std::stringstream ss; ss << comm->getSize() << "ProcLinear" << comm->getRank();
      data_file = fopen(ss.str().c_str(),"r");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file == NULL, MueLu::Exceptions::RuntimeError,"Problem opening file " << ss.str());
      // loop over all local nodes
      fscanf(data_file,"%d",&(dofGlobals[0]));
      fclose(data_file);
      for(decltype(nDofs) i = 0; i < nDofs; i++) dofGlobals[i] = i + dofGlobals[0];
      for(decltype(nDofs) i = 0; i < nDofs; i++) dofGlobals[i] = dofGlobals[i] - 1;
      dofLinearMap = MapFactory::Build(lib,nGlobalDof,dofGlobals(),0,comm);
    }

    Teuchos::RCP<Import> dofImporter = ImportFactory::Build(dofLinearMap, DistributedMap);

    RCP<MultiVector> RHS  = MultiVectorFactory::Build(DistributedMap,1);
    RCP<MultiVector> LHS  = MultiVectorFactory::Build(DistributedMap,1);

    temp = IO::ReadMultiVector ("rhs.mm", dofLinearMap);
    RHS->doImport(*temp, *dofImporter, Xpetra::INSERT);
    LHS->putScalar(Teuchos::ScalarTraits<Scalar>::zero());

    // MueLu part

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

    std::cout << paramList << std::endl;

    const std::string userName = "user data";
    Teuchos::ParameterList& userParamList = paramList.sublist(userName);
    userParamList.set("multivector Coordinates",coordinates);
    userParamList.set("ArrayRCP<LO> DofPresent", dofPresent);

    RCP<Hierarchy> H;
    H = MueLu::CreateXpetraPreconditioner(DistributedMatrix, paramList, paramList /*coordinates*/);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";  clp.setOption("node", &node, "node type (serial | openmp | cuda)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      if (node == "") {
        typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#  else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#  endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_HAVE_SERIAL
        typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_HAVE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_HAVE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
#else
      throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

