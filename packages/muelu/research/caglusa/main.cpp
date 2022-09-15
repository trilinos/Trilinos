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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_RowMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_BlockedMap_def.hpp>

#include <Tpetra_BlockedMatrix_decl.hpp>
#include <Tpetra_BlockedMatrix_def.hpp>

#include <Tpetra_HierarchicalOperator_decl.hpp>
#include <Tpetra_HierarchicalOperator_def.hpp>

#include <Xpetra_TpetraBlockedMap.hpp>
#include <Xpetra_TpetraBlockedMatrix.hpp>

#include <Xpetra_HierarchicalOperator_decl.hpp>
#include <Xpetra_HierarchicalOperator_def.hpp>


#include <Xpetra_IO.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include "MueLu_Exceptions.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct IOhelpers {

  static
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Read(const std::string&   filename,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
       RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap        = Teuchos::null,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rangeMap         = Teuchos::null,
       const bool           callFillComplete = true,
       const bool           binary           = false,
       const bool           readLocal        = false) {
    using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    if (readLocal)
      A = IO::ReadLocal(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    else
      A = IO::Read(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    return A;
  }

  static
  Teuchos::RCP<Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Read(Teuchos::ParameterList& hierarchicalParams,
       RCP< const Teuchos::Comm<int> >& comm) {
    using HOp = Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using blocked_matrix_type = typename HOp::blocked_matrix_type;
    using blocked_map_type = typename blocked_matrix_type::blocked_map_type;
    using matrix_type = typename HOp::matrix_type;
    using map_type = typename HOp::map_type;
    using lo_vec_type = typename blocked_map_type::lo_vec_type;

    auto  lib = Xpetra::UseTpetra;
    RCP<HOp>                       op;
    RCP<const map_type>            map, near_colmap, clusterCoeffMap, ghosted_clusterCoeffMap, clusterMap, ghosted_clusterMap;
    RCP<matrix_type>               nearField, basisMatrix, kernelApproximations, kernelBlockGraph;

    std::vector<RCP<blocked_matrix_type> > transferMatrices;
    RCP<lo_vec_type>               clusterSizes;
    RCP<blocked_map_type>          blockedClusterMap, ghosted_blockedClusterMap;
    RCP<blocked_matrix_type>       blockKernelApproximations;

    const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
    const bool readLocal = hierarchicalParams.get<bool>("read local", false);

    using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    // row, domain and range map of the operator
    map = IO::ReadMap(hierarchicalParams.get<std::string>("map"), lib, comm, readBinary);
    // colmap of near field
    near_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("near colmap"), lib, comm, readBinary);
    // 1-to-1 map for the cluster coefficients
    clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("coefficient map"), lib, comm, readBinary);
    // overlapping map for the cluster coefficients
    ghosted_clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted coefficient map"), lib, comm, readBinary);
    // 1-to-1 map for the clusters
    clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("cluster map"), lib, comm, readBinary);
    // overlapping map for the clusters
    ghosted_clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted cluster map"), lib, comm, readBinary);

    // blocked cluster map
    clusterSizes = Xpetra::IO<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierarchicalParams.get<std::string>("gid_cluster_to_gid_coeff"), clusterMap)->getVectorNonConst(0);
    blockedClusterMap = rcp(new blocked_map_type(clusterCoeffMap, clusterSizes));

    // near field interactions
    nearField = Read(hierarchicalParams.get<std::string>("near field matrix"), map, near_colmap, map, map, true, readBinary, readLocal);

    // far field basis expansion coefficients
    basisMatrix = IOhelpers::Read(hierarchicalParams.get<std::string>("basis expansion coefficient matrix"), map, clusterCoeffMap, clusterCoeffMap, map, true, readBinary, readLocal);

    // far field interactions
    kernelApproximations = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix"), clusterCoeffMap, ghosted_clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
    // block graph of far field interactions
    kernelBlockGraph = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix")+".block", clusterMap, ghosted_clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);

    {
      auto import = kernelBlockGraph->getCrsGraph()->getImporter();
      RCP<lo_vec_type> ghosted_clusterSizes = Xpetra::VectorFactory<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::Build(ghosted_clusterMap);
      ghosted_clusterSizes->doImport(*clusterSizes, *import, Xpetra::INSERT);
      ghosted_blockedClusterMap = rcp(new blocked_map_type(ghosted_clusterCoeffMap, ghosted_clusterSizes));
    }

    blockKernelApproximations = rcp(new blocked_matrix_type(kernelApproximations, kernelBlockGraph, blockedClusterMap, ghosted_blockedClusterMap));

    // Transfer matrices
    auto transfersList = hierarchicalParams.sublist("shift coefficient matrices");
    for (int i = 0; i < transfersList.numParams(); i++) {
      std::string filename = transfersList.get<std::string>(std::to_string(i));
      auto transferPoint = IOhelpers::Read(filename, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
      auto transferBlock = IOhelpers::Read(filename+".block", clusterMap, clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);
      auto transfer = rcp(new blocked_matrix_type(transferPoint, transferBlock, blockedClusterMap));
      transferMatrices.push_back(transfer);
    }

    op = rcp(new HOp(nearField, blockKernelApproximations, basisMatrix, transferMatrices));

    return op;
  }

};


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  #include "MueLu_UseShortNames.hpp"

  std::string xmlHierarchical  = "hierarchical-1d-mm-global.xml";
  std::string xmlMueLu        = "muelu.xml";
  std::string xmlAuxHierarchy = "auxiliary.xml";
  clp.setOption("xml",    &xmlHierarchical);
  clp.setOption("xmlMueLu", &xmlMueLu);
  clp.setOption("xmlAux", &xmlAuxHierarchy);
  bool printTimings  = true; clp.setOption("timings", "notimings", &printTimings,  "print timings to screen");
  bool doTests       = true; clp.setOption("tests",   "notests",   &doTests,       "Test operator using known LHS & RHS.");
  bool doUnPrecSolve = true; clp.setOption("unPrec",  "noUnPrec",  &doUnPrecSolve, "Solve unpreconditioned");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("Hierarchical Driver"));
  Teuchos::RCP<Teuchos::FancyOStream> verbose_out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
  verbose_out->setShowProcRank(true);
  stacked_timer->setVerboseOstream(verbose_out);
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  using HOp = Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using blocked_matrix_type = typename HOp::blocked_matrix_type;
  using blocked_map_type = typename blocked_matrix_type::blocked_map_type;
  using matrix_type = typename HOp::matrix_type;
  using map_type = typename HOp::map_type;
  using mv_type = typename HOp::mv_type;
  using lo_vec_type = typename blocked_map_type::lo_vec_type;
  using coord_mv = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>;
  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using IOhelpers = IOhelpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);
  bool success = true;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const MagnitudeType tol = 10000*Teuchos::ScalarTraits<MagnitudeType>::eps();

  Teuchos::ParameterList         hierarchicalParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlHierarchical, Teuchos::Ptr<Teuchos::ParameterList>(&hierarchicalParams), *comm);

  RCP<HOp> op;
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Read hierarchical matrix")));
    op = IOhelpers::Read(hierarchicalParams, comm);
  }

  out << "Compression: " << op->getCompression() << " of dense matrix."<< std::endl;

  RCP<const map_type> map = op->getDomainMap();
  RCP<matrix_type>    auxOp;
  RCP<mv_type>        X_ex, RHS, X;
  RCP<coord_mv>       coords;
  {
    // Read in auxiliary stuff

    const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
    const bool readLocal = hierarchicalParams.get<bool>("read local", false);

    // Auxiliary matrix used for multigrid construction
    const std::string auxOpStr = hierarchicalParams.get<std::string>("auxiliary operator");
    if (auxOpStr == "near")
      auxOp = op->nearFieldMatrix();
    else {
      // colmap of auxiliary operator
      RCP<const map_type> aux_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("aux colmap"), lib, comm, readBinary);

      auxOp  = IOhelpers::Read(auxOpStr, map, aux_colmap, map, map, true, readBinary, readLocal);
    }

    X_ex = IO::ReadMultiVector(hierarchicalParams.get<std::string>("exact solution"), map);
    RHS  = IO::ReadMultiVector(hierarchicalParams.get<std::string>("right-hand side"), map);
    X    = MultiVectorFactory::Build(map, 1);

    coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierarchicalParams.get<std::string>("coordinates"), map);
  }

  if (doTests) {
    // Some simple apply tests
    Scalar opX_exRHS, MopX_exRHS, MopTX_exRHS;
    {
      op->apply(*X_ex, *X);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
      X->update(one, *RHS, -one);
      opX_exRHS = X->getVector(0)->norm2();
      out << "|op*X_ex - RHS| = " << opX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff.mtx", *X);
    }

    {
      op->apply(*X_ex, *X, Teuchos::NO_TRANS, -one);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
      X->update(one, *RHS, one);
      MopX_exRHS = X->getVector(0)->norm2();
      out << "|(-op)*X_ex + RHS| = " << MopX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
    }

    {
      op->apply(*X_ex, *X, Teuchos::TRANS, -one);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
      X->update(one, *RHS, one);
      MopTX_exRHS = X->getVector(0)->norm2();
      out << "|(-op^T)*X_ex + RHS| = " << MopTX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
    }

    TEUCHOS_ASSERT(opX_exRHS < tol);
    TEUCHOS_ASSERT(MopX_exRHS < tol);
    TEUCHOS_ASSERT(MopTX_exRHS < tol);
  }

#ifdef HAVE_MUELU_BELOS
  if (doUnPrecSolve) {
    // Solve linear system using unpreconditioned Krylov method
    out << "\n*********************************************************\n";
    out << "Unpreconditioned Krylov method\n";
    out << "*********************************************************\n\n";

    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Unpreconditioned solve")));

    using MV = typename HOp::mv_type;
    using OP = Belos::OperatorT<MV>;

    X->putScalar(zero);
    RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

    std::string belosType = "Pseudoblock CG";
    RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
    belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency",      1);
    belosList->set("Output Style",          Belos::Brief);

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType ret = solver->solve();
    int numIts = solver->getNumIters();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << numIts << std::endl;

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *X_ex, -one);
    out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl << std::endl;

    success &= (ret == Belos::Converged);

  }
#endif // HAVE_MUELU_BELOS

  {
    // Solve linear system using a AMG preconditioned Krylov method

    RCP<Hierarchy> auxH, H;

    {
      ////////////////////////////////////////////////////////////////
      // Build the auxiliary hierarchy
      out << "\n*********************************************************\n";
      out << "Building the auxiliary hierarchy\n";
      out << "*********************************************************\n\n";

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Construct auxiliary hierarchy")));

      Teuchos::ParameterList auxParams;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlAuxHierarchy, Teuchos::Ptr<Teuchos::ParameterList>(&auxParams), *comm);
      auxParams.sublist("user data").set("Coordinates", coords);
      // TEUCHOS_ASSERT_EQUALITY(auxParams.get("multigrid algorithm", "unsmoothed"), "unsmoothed");

      auxH = MueLu::CreateXpetraPreconditioner(auxOp, auxParams);
    }

    {
      ////////////////////////////////////////////////////////////////
      // Construct the main hierarchy
      out << "\n*********************************************************\n";
      out << "Building the main hierarchy\n";
      out << "*********************************************************\n\n";

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Construct hierarchy")));

      Teuchos::ParameterList params;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlMueLu, Teuchos::Ptr<Teuchos::ParameterList>(&params), *comm);
      params.set("coarse: max size", 1);
      params.set("max levels", auxH->GetNumLevels());

      op->describe(out, Teuchos::VERB_EXTREME);

      H = rcp(new Hierarchy());
      RCP<Level> lvl = H->GetLevel(0);
      lvl->Set("A", rcp_dynamic_cast<Operator>(op));
      lvl->Set("Coordinates", coords);
      for(int lvlNo = 1; lvlNo < auxH->GetNumLevels(); lvlNo++) {
        H->AddNewLevel();
        RCP<Level> auxLvl = auxH->GetLevel(lvlNo);
        // auto mgr = auxLvl->GetFactoryManager();
        // auxLvl->print(std::cout, MueLu::Debug);
        RCP<Level> fineLvl = H->GetLevel(lvlNo-1);
        lvl = H->GetLevel(lvlNo);
        auto P = auxLvl->Get<RCP<Matrix> >("P");
        auto fineA = rcp_dynamic_cast<HOp>(fineLvl->Get<RCP<Operator> >("A"));
        lvl->Set("P", P);
        params.sublist("level "+std::to_string(lvlNo)).set("P", P);

        auto coarseA = fineA->restrict(P);
        coarseA->describe(out, Teuchos::VERB_EXTREME);
        if (lvlNo+1 == auxH->GetNumLevels())
          lvl->Set("A", coarseA->toMatrix());
        else
          lvl->Set("A", rcp_dynamic_cast<Operator>(coarseA));
      }

      RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(params,op->getDomainMap()->getComm()));
      H->setlib(op->getDomainMap()->lib());
      H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());
      mueLuFactory->SetupHierarchy(*H);
      H->IsPreconditioner(true);
    }


#ifdef HAVE_MUELU_BELOS
    {
      ////////////////////////////////////////////////////////////////
      // Set up the Krylov solver

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Preconditioned solve")));

      using MV = typename HOp::mv_type;
      using OP = Belos::OperatorT<MV>;

      X->putScalar(zero);
      RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
      RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H));
      RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

      std::string belosType = "Pseudoblock CG";
      RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
      belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
      belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
      belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList->set("Output Frequency",      1);
      belosList->set("Output Style",          Belos::Brief);

      belosProblem->setRightPrec(belosPrec);

      bool set = belosProblem->setProblem();
      if (set == false) {
        throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
      }

      // Create an iterative solver manager
      Belos::SolverFactory<Scalar, MV, OP> solverFactory;
      RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
      solver->setProblem(belosProblem);

      // Perform solve
      Belos::ReturnType ret = solver->solve();
      int numIts = solver->getNumIters();

      // Get the number of iterations for this solve.
      out << "Number of iterations performed for this solve: " << numIts << std::endl;

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
      X->update(one, *X_ex, -one);
      out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl;

      success &= (ret == Belos::Converged);
    }

    stacked_timer->stop("Hierarchical Driver");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    if (printTimings)
      stacked_timer->report(out, comm, options);

#endif // HAVE_MUELU_BELOS
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
