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

#include <Xpetra_HierarchicalOperator_decl.hpp>
#include <Xpetra_HierarchicalOperator_def.hpp>

#include <Xpetra_IO.hpp>
#include <MueLu_IOhelpers.hpp>
#include <auxiliaryOperators.hpp>

#include <MueLu.hpp>
#include "MueLu_Exceptions.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#endif

#define MUELU_HIERARCHICAL_DEBUG

using Teuchos::RCP;
using Teuchos::rcp;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  std::string xmlHierarchical = "1d-binary/hierarchical.xml";
  clp.setOption("xmlHierarchical", &xmlHierarchical, "XML describing the hierarchical operator");
  std::string xmlProblem = "1d-binary/problem.xml";
  clp.setOption("xmlProblem", &xmlProblem, "XML describing the problem");
  std::string xmlBelos = "belos.xml";
  clp.setOption("xmlBelos", &xmlBelos, "XML with Belos parameters");
  std::string xmlMueLu = "muelu.xml";
  clp.setOption("xmlMueLu", &xmlMueLu, "XML with MueLu parameters");
  std::string xmlAuxHierarchy = "auxiliary.xml";
  clp.setOption("xmlAux", &xmlAuxHierarchy, "XML with MueLu parameters for the auxiliary hierarchy");
  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  bool doTests = true;
  clp.setOption("tests", "notests", &doTests, "Test operator using known LHS & RHS.");
  bool doUnPrecSolve = true;
  clp.setOption("unPrec", "noUnPrec", &doUnPrecSolve, "Solve unpreconditioned");
  bool doPrecSolve = true;
  clp.setOption("prec", "noPrec", &doPrecSolve, "Solve preconditioned with AMG");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("Hierarchical Driver"));
  Teuchos::RCP<Teuchos::FancyOStream> verbose_out   = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
  verbose_out->setShowProcRank(true);
  stacked_timer->setVerboseOstream(verbose_out);
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  using HOp                 = Xpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using op_type             = Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using blocked_matrix_type = typename HOp::blocked_matrix_type;
  using blocked_map_type    = typename blocked_matrix_type::blocked_map_type;
  using matrix_type         = typename HOp::matrix_type;
  using map_type            = typename HOp::map_type;
  using mv_type             = typename HOp::mv_type;
  using lo_vec_type         = typename blocked_map_type::lo_vec_type;
  using coord_mv            = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>;
  using MagnitudeType       = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using IO                  = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using IOhelpers           = MueLu::IOhelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream &out       = *fancy;
  out.setOutputToRootOnly(0);
  bool success            = true;
  const Scalar one        = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero       = Teuchos::ScalarTraits<Scalar>::zero();
  const MagnitudeType tol = 100000 * Teuchos::ScalarTraits<MagnitudeType>::eps();

  RCP<op_type> op;
  RCP<HOp> hop;
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Read hierarchical matrix")));

    op  = IOhelpers::Read(xmlHierarchical, comm);
    hop = Teuchos::rcp_dynamic_cast<HOp>(op);
  }

  if (!hop.is_null())
    out << "Compression: " << hop->getCompression() << " of dense matrix." << std::endl;

  Teuchos::ParameterList problemParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlProblem, Teuchos::Ptr<Teuchos::ParameterList>(&problemParams), *comm);

  RCP<const map_type> map = op->getDomainMap();
  RCP<matrix_type> auxOp, mass;
  RCP<mv_type> X_ex, RHS, X;
  RCP<coord_mv> coords;
  {
    // Read in auxiliary stuff

    // coordinates
    coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(problemParams.get<std::string>("coordinates"), map);

    // Auxiliary matrix used for multigrid construction
    {
      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Construct auxiliary operator")));

      auxOp = MueLu::constructAuxiliaryOperator(op, problemParams);
    }

    // Mass matrix for L2 error computation
    {
      const bool readBinary = problemParams.get<bool>("read binary", false);
      const bool readLocal  = problemParams.get<bool>("read local", false);
      // colmap of auxiliary operator
      RCP<const map_type> aux_colmap = IO::ReadMap(problemParams.get<std::string>("aux colmap"), lib, comm, readBinary);

      mass = IOhelpers::Read(problemParams.get<std::string>("mass matrix"), map, aux_colmap, map, map, true, readBinary, readLocal);
    }

    // known pair of LHS, RHS
    X_ex = IO::ReadMultiVector(problemParams.get<std::string>("exact solution"), map);
    RHS  = IO::ReadMultiVector(problemParams.get<std::string>("right-hand side"), map);
    // solution vector
    X = MultiVectorFactory::Build(map, 1);
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
  Teuchos::ParameterList belosParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlBelos, Teuchos::Ptr<Teuchos::ParameterList>(&belosParams), *comm);
  if (doUnPrecSolve) {
    // Solve linear system using unpreconditioned Krylov method
    out << "\n*********************************************************\n";
    out << "Unpreconditioned Krylov method\n";
    out << "*********************************************************\n\n";

    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Unpreconditioned solve")));

    using MV = typename HOp::mv_type;
    using OP = Belos::OperatorT<MV>;

    X->putScalar(zero);
    RCP<OP> belosOp                                         = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

    std::string belosType = "Pseudoblock CG";
    auto belosSolverList  = rcpFromRef(belosParams.sublist(belosType));

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP<Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosSolverList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType ret = solver->solve();
    int numIts            = solver->getNumIters();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << numIts << std::endl;

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *X_ex, -one);
    out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl
        << std::endl;

    success &= (ret == Belos::Converged);
  }
#endif  // HAVE_MUELU_BELOS

  if (doPrecSolve) {
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
      auxParams.set("hierarchy label", "Auxiliary");
      auxParams.sublist("user data").set("Coordinates", coords);
      // No rebalancing yet.
      auxParams.set("coarse: max size", std::max(auxParams.get("coarse: max size", 2 * comm->getSize()),
                                                 2 * comm->getSize()));

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
      params.sublist("user data").set("Coordinates", coords);

      H = MueLu::constructHierarchyFromAuxiliary(Teuchos::rcp_dynamic_cast<Operator>(op, true), auxH, params, out);
    }

#ifdef HAVE_MUELU_BELOS
    {
      ////////////////////////////////////////////////////////////////
      // Set up the Krylov solver

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Preconditioned solve")));

      using MV = typename HOp::mv_type;
      using OP = Belos::OperatorT<MV>;

      X->putScalar(zero);
      RCP<OP> belosOp                                         = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
      RCP<OP> belosPrec                                       = rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));
      RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

      std::string belosType = "Pseudoblock CG";
      auto belosSolverList  = rcpFromRef(belosParams.sublist(belosType));

      belosProblem->setRightPrec(belosPrec);

      bool set = belosProblem->setProblem();
      if (set == false) {
        throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
      }

      // Create an iterative solver manager
      Belos::SolverFactory<Scalar, MV, OP> solverFactory;
      RCP<Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosSolverList);
      solver->setProblem(belosProblem);

      // Perform solve
      Belos::ReturnType ret = solver->solve();
      int numIts            = solver->getNumIters();

      // Get the number of iterations for this solve.
      out << "Number of iterations performed for this solve: " << numIts << std::endl;

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
      X->update(one, *X_ex, -one);
      out << "|X-X_ex|_2 = " << X->getVector(0)->norm2() << std::endl;

      auto massDiff = MultiVectorFactory::Build(map, 1);
      mass->apply(*X, *massDiff);
      out << "|x-x_ex|_L2 = " << Teuchos::ScalarTraits<Scalar>::squareroot(X->getVector(0)->dot(*massDiff->getVector(0))) << std::endl;

      success &= (ret == Belos::Converged);
    }
#endif  // HAVE_MUELU_BELOS
  }

  stacked_timer->stop("Hierarchical Driver");
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = options.output_histogram = options.output_minmax = true;
  if (printTimings)
    stacked_timer->report(out, comm, options);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
