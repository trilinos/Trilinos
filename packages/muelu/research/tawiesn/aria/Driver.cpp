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

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp
#include <MueLu_VisualizationHelpers.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ExportVTK : public MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  ExportVTK(){};

 public:
  void writeFile(std::ofstream& fout, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coordinates, Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& sol) {
    using namespace std;
    typedef MueLu::VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> VH;
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > nodeMap = coordinates->getMap();
    std::vector<LocalOrdinal> vertices;
    std::vector<LocalOrdinal> geomSize;
    LocalOrdinal numFineNodes = Teuchos::as<LocalOrdinal>(coordinates->getLocalLength());

    vertices.reserve(numFineNodes);
    geomSize.reserve(numFineNodes);
    for (LocalOrdinal i = 0; i < numFineNodes; i++) {
      vertices.push_back(i);
      geomSize.push_back(1);
    }

    Teuchos::ArrayRCP<const double> xCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(0));
    Teuchos::ArrayRCP<const double> yCoords = Teuchos::arcp_reinterpret_cast<const double>(coordinates->getData(1));
    Teuchos::ArrayRCP<const double> zCoords = Teuchos::null;
    if (coordinates->getNumVectors() == 3) {
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
    for (int i = 0; i < int(uniqueFine.size()); i++) {
      fout << myRank << " ";
      if (i % 20 == 19)
        fout << std::endl
             << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    // solution vector
    fout << "        <DataArray type=\"Float64\" Name=\"Solution\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    fout << indent;
    for (int i = 0; i < int(uniqueFine.size()); i++) {
      size_t numVec = coordinates->getNumVectors();
      for (int d = 0; d < numVec; d++)
        fout << solData[i * numVec + d] << " ";
      if (i % 3 == 0)
        fout << std::endl
             << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "      </PointData>" << std::endl;  // that is annoying

    this->writeFileVTKCoordinates(fout, uniqueFine, xCoords, yCoords, zCoords, coordinates->getNumVectors());
    this->writeFileVTKCells(fout, uniqueFine, vertices, geomSize);
    this->writeFileVTKClosing(fout);
    fout.close();
  };
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;

  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    SC zero = Teuchos::ScalarTraits<SC>::zero();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out       = *fancy;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================

    // GO nx = 100, ny = 100, nz = 100;
    // Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of Xpetra

    std::string xmlFileName = "driver.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    // =========================================================================
    // Read in problem
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;

    GlobalOrdinal nGlobalNodes = 0;
    if (comm->getRank() == 0) {
      std::ifstream data_file;
      data_file.open("xxx.mm");
      TEUCHOS_TEST_FOR_EXCEPTION(data_file.good() == false, MueLu::Exceptions::RuntimeError, "Problem opening file xxx.mm");
      std::string line;
      getline(data_file, line);
      data_file >> nGlobalNodes;
      data_file.close();
    }
    Teuchos::broadcast(*comm, 0, 1, &nGlobalNodes);

    out << "Found " << nGlobalNodes << " nodes " << std::endl;

    // each processor reads in its dofPresent array
    LocalOrdinal nNodes = 0;
    LocalOrdinal nDofs  = 0;
    int maxDofPerNode   = -1;
    Teuchos::ArrayRCP<LocalOrdinal> dofPresent;
    {
      std::stringstream ss;
      ss << comm->getSize() << "dofPresent" << comm->getRank();

      try {
        std::ifstream data_file(ss.str());
        // read in first line containing number of local nodes and maxDofPerNode
        data_file >> nNodes >> maxDofPerNode;
        dofPresent = Teuchos::ArrayRCP<LocalOrdinal>(nNodes * maxDofPerNode, 0);
        // loop over all local nodes
        for (LocalOrdinal i = 0; i < nNodes; i++) {
          for (int j = 0; j < maxDofPerNode; j++) {
            int tmp = -1;
            data_file >> tmp;
            if (tmp == 1) {
              dofPresent[i * maxDofPerNode + j] = 1;
              nDofs++;
            } else {
              dofPresent[i * maxDofPerNode + j] = 0;
            }
          }
        }
      } catch (const std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Problem opening/reading file " << ss.str());
      }
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
      std::stringstream ss;
      ss << comm->getSize() << "proc" << comm->getRank();
      try {
        std::ifstream data_file(ss.str());

        // read in first line containing number of local nodes and maxDofPerNode
        LocalOrdinal tmpDofs = 0;
        data_file >> tmpDofs;
        TEUCHOS_TEST_FOR_EXCEPTION(tmpDofs != nDofs, MueLu::Exceptions::RuntimeError, "Number of gid entries in map file is " << tmpDofs << " but should be " << nDofs);

        dofGlobals = Teuchos::Array<GlobalOrdinal>(nDofs);

        // loop over all local nodes
        for (GlobalOrdinal i = 0; i < nDofs; i++) {
          int data;
          data_file >> data;
          dofGlobals[i] = Teuchos::as<GlobalOrdinal>(data);
        }
      } catch (const std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Problem opening/reading file " << ss.str());
      }

      nodalGlobals = Teuchos::Array<GlobalOrdinal>(nNodes);

      GlobalOrdinal count = nDofs - 1;
      for (GlobalOrdinal i = nNodes - 1; i >= 0; i--) {
        for (int j = maxDofPerNode - 1; j >= 0; j--) {
          if (dofPresent[i * maxDofPerNode + j] == 1)
            nodalGlobals[i] = dofGlobals[count--];
        }
      }
    }

    Teuchos::RCP<Map> LapMap = MapFactory::Build(lib, nGlobalNodes, nodalGlobals(), 0, comm);

    // LapMap->describe(out, Teuchos::VERB_EXTREME);

    //

    {
      std::stringstream ss;
      ss << comm->getSize() << "ProcLinear" << comm->getRank();
      try {
        std::ifstream data_file(ss.str());

        // loop over all local nodes
        for (GlobalOrdinal i = 0; i < nNodes; i++) {
          int data;
          data_file >> data;
          nodalGlobals[i] = Teuchos::as<GlobalOrdinal>(data);
        }
      } catch (const std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Problem opening/reading file " << ss.str());
      }
      for (GlobalOrdinal i = 0; i < nNodes; i++) {
        nodalGlobals[i] = nodalGlobals[i] - 1;
      }
    }

    Teuchos::RCP<Map> LinearMap = MapFactory::Build(lib, nGlobalNodes, nodalGlobals(), 0, comm);
    // LinearMap->describe(out, Teuchos::VERB_EXTREME);LapMap->describe(out, Teuchos::VERB_EXTREME);

    // LapMap->describe(out, Teuchos::VERB_EXTREME);LapMap->describe(out, Teuchos::VERB_EXTREME);

    Teuchos::RCP<Import> Importer = ImportFactory::Build(LinearMap, LapMap);

    Teuchos::RCP<MultiVector> xpetraXXX = MultiVectorFactory::Build(LapMap, 1);
    Teuchos::RCP<MultiVector> xpetraYYY = MultiVectorFactory::Build(LapMap, 1);

    RCP<MultiVector> temp = IO::ReadMultiVector("xxx.mm", LinearMap);
    xpetraXXX->doImport(*temp, *Importer, Xpetra::INSERT);
    temp = IO::ReadMultiVector("yyy.mm", LinearMap);
    xpetraYYY->doImport(*temp, *Importer, Xpetra::INSERT);

    Teuchos::RCP<MultiVector> coordinates = MultiVectorFactory::Build(LapMap, 2);
    Teuchos::ArrayRCP<const Scalar> srcX  = xpetraXXX->getData(0);
    Teuchos::ArrayRCP<const Scalar> srcY  = xpetraYYY->getData(0);
    Teuchos::ArrayRCP<Scalar> dataX       = coordinates->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> dataY       = coordinates->getDataNonConst(1);
    for (decltype(coordinates->getLocalLength()) i = 0; i < coordinates->getLocalLength(); i++) {
      dataX[i] = srcX[i];
      dataY[i] = srcY[i];
    }

    // read in matrix
    Teuchos::RCP<const Map> DistributedMap = Teuchos::null;
    Teuchos::RCP<Matrix> DistributedMatrix = Teuchos::null;

    // read in matrix to determine number of rows
    Teuchos::RCP<Matrix> SerialMatrix = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("theMatrix.m", lib, comm);
    GlobalOrdinal NumRows             = SerialMatrix->getRowMap()->getGlobalNumElements();

    // re-read in matrix and distribute it using the user-given distribution over processors
    DistributedMap    = MapFactory::Build(lib, NumRows, dofGlobals(), 0, comm);
    DistributedMatrix = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("theMatrix.m",
                                                                                    DistributedMap,
                                                                                    Teuchos::null,
                                                                                    DistributedMap,
                                                                                    DistributedMap,
                                                                                    true,    //           callFillComplete = true,
                                                                                    false,   //           binary           = false,
                                                                                    false,   //           tolerant         = false,
                                                                                    false);  //           debug            = false)

    // read in global vectors (e.g. rhs)
    GlobalOrdinal nGlobalDof = 0;
    GlobalOrdinal nLocalDofs = Teuchos::as<GlobalOrdinal>(nDofs);

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, comm->getSize(), &nLocalDofs, &nGlobalDof);

    Teuchos::RCP<const Map> dofLinearMap = Teuchos::null;
    {
      std::stringstream ss;
      ss << comm->getSize() << "ProcLinear" << comm->getRank();
      try {
        std::ifstream data_file(ss.str());
        // loop over all local nodes
        data_file >> dofGlobals[0];
      } catch (const std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Problem opening/reading file " << ss.str());
      }
      for (decltype(nDofs) i = 0; i < nDofs; i++) dofGlobals[i] = i + dofGlobals[0];
      for (decltype(nDofs) i = 0; i < nDofs; i++) dofGlobals[i] = dofGlobals[i] - 1;
      dofLinearMap = MapFactory::Build(lib, nGlobalDof, dofGlobals(), 0, comm);
    }

    Teuchos::RCP<Import> dofImporter = ImportFactory::Build(dofLinearMap, DistributedMap);

    RCP<MultiVector> RHS = MultiVectorFactory::Build(DistributedMap, 1);
    RCP<MultiVector> LHS = MultiVectorFactory::Build(DistributedMap, 1);

    temp = IO::ReadMultiVector("rhs.mm", dofLinearMap);
    RHS->doImport(*temp, *dofImporter, Xpetra::INSERT);
    LHS->putScalar(zero);

    // MueLu part

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

    const std::string userName            = "user data";
    Teuchos::ParameterList& userParamList = paramList.sublist(userName);
    userParamList.set("multivector Coordinates", coordinates);
    userParamList.set("ArrayRCP<LO> DofPresent", dofPresent);

    RCP<Hierarchy> H;
    H = MueLu::CreateXpetraPreconditioner(DistributedMatrix, paramList);

#ifdef HAVE_MUELU_BELOS
    {
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Belos Solve")));
      // Operator and Multivector type that will be used with Belos
      typedef MultiVector MV;
      typedef Belos::OperatorT<MV> OP;
      H->IsPreconditioner(true);

      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(DistributedMatrix));  // Turns a Xpetra::Matrix object into a Belos operator
      Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));                   // Turns a MueLu::Hierarchy object into a Belos operator
      // Construct a Belos LinearProblem object
      RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, LHS, RHS));
      belosProblem->setRightPrec(belosPrec);
      bool set = belosProblem->setProblem();
      if (set == false) {
        out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos parameter list
      Teuchos::ParameterList belosList;
      belosList.set("Maximum Iterations", 300);      // Maximum number of iterations allowed
      belosList.set("Convergence Tolerance", 1e-8);  // Relative convergence tolerance requested
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList.set("Output Frequency", 1);
      belosList.set("Output Style", Belos::Brief);

      // Create an iterative solver manager
      RCP<Belos::SolverManager<SC, MV, OP> > solver;
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

      // Perform solve
      Belos::ReturnType ret = Belos::Unconverged;
      try {
        ret = solver->solve();

        // Get the number of iterations for this solve.
        out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      } catch (...) {
        out << std::endl
            << "ERROR:  Belos threw an error! " << std::endl;
      }
      // Check convergence
      if (ret != Belos::Converged)
        out << std::endl
            << "ERROR:  Belos did not converge or a warning occured! " << std::endl;
      else
        out << std::endl
            << "SUCCESS:  Belos converged!" << std::endl;
    }
#endif

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";
    clp.setOption("node", &node, "node type (serial | openmp | cuda | hip)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double, int, int, Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
      if (node == "") {
        typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_ENABLE_SERIAL
        typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_ENABLE_OPENMP
        typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_ENABLE_CUDA
        typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else if (node == "hip") {
#ifdef KOKKOS_ENABLE_HIP
        typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("HIP node type is disabled");
#endif
      } else if (node == "sycl") {
#ifdef KOKKOS_ENABLE_SYCL
        typedef Tpetra::KokkosCompatKokkosSYCLWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("SYCL node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
