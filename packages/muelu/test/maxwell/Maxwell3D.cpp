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
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_IO.hpp>

// MueLu
#include <MueLu_RefMaxwell.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_TestHelpers_Common.hpp>
#include <MueLu_Exceptions.hpp>

#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#endif

// Belos
#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
//#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

// Stratimikos
#ifdef HAVE_MUELU_STRATIMIKOS
// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>
#endif

// Main wrappers struct
// Because C++ doesn't support partial template specialization of functions.
// By default, do not try to run Stratimikos, since that only works for Scalar=double.
template<typename Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
struct MainWrappers {
  static int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]);
};


// Partial template specialization on SC=double
// This code branch gives the option to run with Stratimikos.
template<class LocalOrdinal, class GlobalOrdinal, class Node>
struct MainWrappers<double,LocalOrdinal,GlobalOrdinal,Node> {
  static int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]);
};


// By default, do not try to run Stratimikos, since that only works for Scalar=double.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int MainWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    bool        printTimings      = true;              clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
    std::string timingsFormat     = "table-fixed";     clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
    double scaling                = 1.0;               clp.setOption("scaling",               &scaling,           "scale mass term");
    std::string solverName = "Belos";                  clp.setOption("solverName",            &solverName,        "Name of iterative linear solver "
                                                                                                                  "to use for solving the linear system. "
                                                                                                                  "(\"Belos\")");
    std::string belosSolverType   = "Block CG";        clp.setOption("belosSolverType",       &belosSolverType,   "Name of the Belos linear solver");
    bool        usePrec           = true;              clp.setOption("usePrec", "noPrec",     &usePrec,           "use RefMaxwell preconditioner");
    std::string xml               = "";                clp.setOption("xml",                   &xml,               "xml file with solver parameters");
    double      tol               = 1e-10;             clp.setOption("tol",                   &tol,               "solver convergence tolerance");
    
    std::string S_file, SM_file, M1_file, M0_file, M0inv_file, D0_file, coords_file, rhs_file="", nullspace_file="";

    if (!TYPE_EQUAL(SC, std::complex<double>)) {
      S_file = "S.mat";
      SM_file = "";
      M1_file = "M1.mat";
      M0_file = "M0.mat";
      M0inv_file = "";
      D0_file = "D0.mat";
    } else {
      S_file = "S_complex.mat";
      SM_file = "";
      M1_file = "M1_complex.mat";
      M0_file = "M0_complex.mat";
      M0inv_file = "";
      D0_file = "D0_complex.mat";
    }
    coords_file = "coords.mat";

    clp.setOption("S", &S_file);
    clp.setOption("SM", &SM_file);
    clp.setOption("M1", &M1_file);
    clp.setOption("M0", &M0_file);
    clp.setOption("M0inv", &M0inv_file);
    clp.setOption("D0", &D0_file);
    clp.setOption("coords", &coords_file);
    clp.setOption("nullspace", &nullspace_file);
    clp.setOption("rhs", &rhs_file);

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    comm->barrier();
    auto globalTimeMonitor = TimeMonitor::getNewTimer("Maxwell: S - Global Time");
    auto tm                = TimeMonitor::getNewTimer("Maxwell: 1 - Read and Build Matrices");

    // Read matrices in from files
    // gradient matrix
    RCP<Matrix> D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    // maps for nodal and edge matrices
    RCP<const Map> node_map = D0_Matrix->getDomainMap();
    RCP<const Map> edge_map = D0_Matrix->getRangeMap();
    // edge mass matrix
    RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
    // build stiffness plus mass matrix (SM_Matrix)
    RCP<Matrix> SM_Matrix;
    if (SM_file == "") {
      // edge stiffness matrix
      RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*S_Matrix,false,(SC)1.0,*M1_Matrix,false,scaling,SM_Matrix,*out);
      SM_Matrix->fillComplete();
    } else {
      SM_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(SM_file, edge_map);
    }
    RCP<Matrix> M0inv_Matrix;
    if (M0inv_file == "") {
      // nodal mass matrix
      RCP<Matrix> M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0_file, node_map);
      // build lumped mass matrix inverse (M0inv_Matrix)
      RCP<Vector> diag = Utilities::GetLumpedMatrixDiagonal(M0_Matrix);
      RCP<CrsMatrixWrap> M0inv_MatrixWrap = Teuchos::rcp(new CrsMatrixWrap(node_map, node_map, 0, Xpetra::StaticProfile));
      RCP<CrsMatrix> M0inv_CrsMatrix = M0inv_MatrixWrap->getCrsMatrix();
      Teuchos::ArrayRCP<size_t> rowPtr;
      Teuchos::ArrayRCP<LO> colInd;
      Teuchos::ArrayRCP<SC> values;
      Teuchos::ArrayRCP<const SC> diags = diag->getData(0);
      size_t nodeNumElements = node_map->getNodeNumElements();
      M0inv_CrsMatrix->allocateAllValues(nodeNumElements, rowPtr, colInd, values);
      SC ONE = (SC)1.0;
      for (size_t i = 0; i < nodeNumElements; i++) {
        rowPtr[i] = i;  colInd[i] = i;  values[i] = ONE / diags[i];
      }
      rowPtr[nodeNumElements] = nodeNumElements;
      M0inv_CrsMatrix->setAllValues(rowPtr, colInd, values);
      M0inv_CrsMatrix->expertStaticFillComplete(node_map, node_map);
      M0inv_Matrix = Teuchos::rcp_dynamic_cast<Matrix>(M0inv_MatrixWrap);
    } else if (M0inv_file == "none") {
      // pass
    } else {
      M0inv_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0inv_file, node_map);
    }
    // coordinates
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> > coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

    RCP<MultiVector> nullspace = Teuchos::null;
    if (nullspace_file != "")
      nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace_file, edge_map);

    // set parameters
    std::string defaultXMLfile;
    if (!TYPE_EQUAL(SC, std::complex<double>) && !TYPE_EQUAL(SC, std::complex<float>))
      defaultXMLfile = "Maxwell.xml";
    else
      defaultXMLfile = "Maxwell_complex.xml";
    Teuchos::ParameterList params;
    Teuchos::updateParametersFromXmlFileAndBroadcast(defaultXMLfile,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);
    if (xml != "")
      Teuchos::updateParametersFromXmlFileAndBroadcast(xml,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);

    // setup LHS, RHS
    RCP<MultiVector> B;
    if (rhs_file == "") {
      B = MultiVectorFactory::Build(edge_map,1);
      RCP<MultiVector> vec = MultiVectorFactory::Build(edge_map,1);
      vec -> putScalar((SC)1.0);
      SM_Matrix->apply(*vec,*B);
    } else
      B = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhs_file, edge_map);

    RCP<MultiVector> X = MultiVectorFactory::Build(edge_map,1);
    X -> putScalar((SC)0.0);

    comm->barrier();
    tm = Teuchos::null;

    if (solverName == "Belos") {
      auto tm2  = TimeMonitor::getNewTimer("Maxwell: 2 - Build Belos solver etc");

      // construct preconditioner
      RCP<MueLu::RefMaxwell<SC,LO,GO,NO> > preconditioner;
      if (usePrec)
        preconditioner
          = rcp( new MueLu::RefMaxwell<SC,LO,GO,NO>(SM_Matrix,D0_Matrix,M0inv_Matrix,
                                                    M1_Matrix,nullspace,coords,params) );


#ifdef HAVE_MUELU_TPETRA
      {
        // A test to make sure we can wrap this guy as a MueLu::TpetraOperator
        RCP<Operator> precOp = Teuchos::rcp_dynamic_cast<Operator>(preconditioner);
        MueLu::TpetraOperator<SC,LO,GO,NO> OpT(precOp);
      }
#endif

      // Belos linear problem
#ifdef HAVE_MUELU_BELOS
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(SM_Matrix)); // Turns a Xpetra::Matrix object into a Belos operator

      RCP<Belos::LinearProblem<SC, MV, OP> > problem = rcp( new Belos::LinearProblem<SC, MV, OP>() );
      problem -> setOperator( belosOp );
      Teuchos::RCP<OP> belosPrecOp;
      if (usePrec) {
        belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner)); // Turns a Xpetra::Matrix object into a Belos operator
        problem -> setRightPrec( belosPrecOp );
      }

      problem -> setProblem( X, B );

      bool set = problem->setProblem();
      if (set == false) {
        *out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos solver
      RCP< Belos::SolverManager<SC, MV, OP> > solver;
      RCP< Belos::SolverFactory<SC, MV,OP> > factory = rcp( new  Belos::SolverFactory<SC,MV,OP>() );
      RCP<Teuchos::ParameterList> belosParams
        = rcp( new Teuchos::ParameterList() );
      belosParams->set("Maximum Iterations", 100);
      belosParams->set("Convergence Tolerance",tol);
      belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosParams->set("Output Frequency",1);
      belosParams->set("Output Style",Belos::Brief);
      solver = factory->create(belosSolverType,belosParams);

      comm->barrier();
      tm2=Teuchos::null;
      auto tm3  = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // set problem and solve
      solver -> setProblem( problem );
      Belos::ReturnType status = solver -> solve();
      int iters = solver -> getNumIters();
      success = (iters<50 && status == Belos::Converged);
      if (success)
        *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        *out << "FAILURE! Belos did not converge fast enough." << std::endl;
      tm3 = Teuchos::null;
    }
    comm->barrier();
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
      if (timingsFormat == "yaml") {
        reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
        reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
      }
      reportParams->set("How to merge timer sets",   "Union");
      reportParams->set("alwaysWriteLocal",          false);
      reportParams->set("writeGlobalStats",          true);
      reportParams->set("writeZeroTimers",           false);
      // FIXME: no "ignoreZeroTimers"

      const std::string filter = "";

      std::ios_base::fmtflags ff(out->flags());
      if (timingsFormat == "table-fixed") *out << std::fixed;
      else * out << std::scientific;
      TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
       *out << std::setiosflags(ff);
    }

    TimeMonitor::clearCounters();
#endif // #ifdef HAVE_MUELU_BELOS
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#else
  return EXIT_SUCCESS;
#endif
} // main


// This code branch gives the option to run with Stratimikos.
template<class LocalOrdinal, class GlobalOrdinal, class Node>
int MainWrappers<double,LocalOrdinal,GlobalOrdinal,Node>::main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  typedef double Scalar;
#include <MueLu_UseShortNames.hpp>

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

// #if defined(HAVE_TPETRA_INST_INT_INT)

#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    bool        printTimings      = true;              clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
    std::string timingsFormat     = "table-fixed";     clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
    double scaling                = 1.0;               clp.setOption("scaling",               &scaling,           "scale mass term");
    std::string solverName        = "Belos";           clp.setOption("solverName",            &solverName,        "Name of iterative linear solver "
                                                                                                                  "to use for solving the linear system. "
                                                                                                                  "(\"Belos\" or \"Stratimikos\")");
    std::string belosSolverType   = "Block CG";        clp.setOption("belosSolverType",       &belosSolverType,   "Name of the Belos linear solver");
    bool        usePrec           = true;              clp.setOption("usePrec", "noPrec",     &usePrec,           "use RefMaxwell preconditioner");
    std::string xml               = "";                clp.setOption("xml",                   &xml,               "xml file with solver parameters");
    double      tol               = 1e-10;             clp.setOption("tol",                   &tol,               "solver convergence tolerance");

    std::string S_file, SM_file, M1_file, M0_file, M0inv_file, D0_file, coords_file, rhs_file="", nullspace_file="";
    if (!TYPE_EQUAL(SC, std::complex<double>)) {
      S_file = "S.mat";
      SM_file = "";
      M1_file = "M1.mat";
      M0_file = "M0.mat";
      M0inv_file = "";
      D0_file = "D0.mat";
    } else {
      S_file = "S_complex.mat";
      SM_file = "";
      M1_file = "M1_complex.mat";
      M0_file = "M0_complex.mat";
      M0inv_file = "";
      D0_file = "D0_complex.mat";
    }
    coords_file = "coords.mat";

    clp.setOption("S", &S_file);
    clp.setOption("SM", &SM_file);
    clp.setOption("M1", &M1_file);
    clp.setOption("M0", &M0_file);
    clp.setOption("M0inv", &M0inv_file);
    clp.setOption("D0", &D0_file);
    clp.setOption("coords", &coords_file);
    clp.setOption("rhs", &rhs_file);
    clp.setOption("nullspace", &nullspace_file);

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    comm->barrier();
    auto globalTimeMonitor = TimeMonitor::getNewTimer("Maxwell: S - Global Time");
    auto tm                = TimeMonitor::getNewTimer("Maxwell: 1 - Read and Build Matrices");

    // Read matrices in from files
    // gradient matrix
    RCP<Matrix> D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    // maps for nodal and edge matrices
    RCP<const Map> node_map = D0_Matrix->getDomainMap();
    RCP<const Map> edge_map = D0_Matrix->getRangeMap();
    // edge mass matrix
    RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
    // build stiffness plus mass matrix (SM_Matrix)
    RCP<Matrix> SM_Matrix;
    if (SM_file == "") {
      // edge stiffness matrix
      RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*S_Matrix,false,(SC)1.0,*M1_Matrix,false,scaling,SM_Matrix,*out);
      SM_Matrix->fillComplete();
    } else {
      SM_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(SM_file, edge_map);
    }
    RCP<Matrix> M0inv_Matrix;
    if (M0inv_file == "") {
      // nodal mass matrix
      RCP<Matrix> M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0_file, node_map);
      // build lumped mass matrix inverse (M0inv_Matrix)
      RCP<Vector> diag = Utilities::GetLumpedMatrixDiagonal(M0_Matrix);
      RCP<CrsMatrixWrap> M0inv_MatrixWrap = Teuchos::rcp(new CrsMatrixWrap(node_map, node_map, 0, Xpetra::StaticProfile));
      RCP<CrsMatrix> M0inv_CrsMatrix = M0inv_MatrixWrap->getCrsMatrix();
      Teuchos::ArrayRCP<size_t> rowPtr;
      Teuchos::ArrayRCP<LO> colInd;
      Teuchos::ArrayRCP<SC> values;
      Teuchos::ArrayRCP<const SC> diags = diag->getData(0);
      size_t nodeNumElements = node_map->getNodeNumElements();
      M0inv_CrsMatrix->allocateAllValues(nodeNumElements, rowPtr, colInd, values);
      SC ONE = (SC)1.0;
      for (size_t i = 0; i < nodeNumElements; i++) {
        rowPtr[i] = i;  colInd[i] = i;  values[i] = ONE / diags[i];
      }
      rowPtr[nodeNumElements] = nodeNumElements;
      M0inv_CrsMatrix->setAllValues(rowPtr, colInd, values);
      M0inv_CrsMatrix->expertStaticFillComplete(node_map, node_map);
      M0inv_Matrix = Teuchos::rcp_dynamic_cast<Matrix>(M0inv_MatrixWrap);
    } else if (M0inv_file == "none") {
      // pass
    } else {
      M0inv_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0inv_file, node_map);
    }
    // coordinates
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

    RCP<MultiVector> nullspace = Teuchos::null;
    if (nullspace_file != "")
      nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace_file, edge_map);

    // set parameters
    std::string defaultXMLfile;
    if (!TYPE_EQUAL(SC, std::complex<double>))
      defaultXMLfile = "Maxwell.xml";
    else
      defaultXMLfile = "Maxwell_complex.xml";
    Teuchos::ParameterList params;
    Teuchos::updateParametersFromXmlFileAndBroadcast(defaultXMLfile,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);
    if (xml != "")
      Teuchos::updateParametersFromXmlFileAndBroadcast(xml,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);

    // setup LHS, RHS
    // setup LHS, RHS
    RCP<MultiVector> B;
    if (rhs_file == "") {
      B = MultiVectorFactory::Build(edge_map,1);
      RCP<MultiVector> vec = MultiVectorFactory::Build(edge_map,1);
      vec -> putScalar((SC)1.0);
      SM_Matrix->apply(*vec,*B);
    } else
      B = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhs_file, edge_map);
    RCP<MultiVector> X = MultiVectorFactory::Build(edge_map,1);
    X -> putScalar((SC)0.0);

    comm->barrier();
    tm = Teuchos::null;

    if (solverName == "Belos") {
      auto tm2 = TimeMonitor::getNewTimer("Maxwell: 2 - Build Belos solver etc");

      // construct preconditioner
      RCP<MueLu::RefMaxwell<SC,LO,GO,NO> > preconditioner;
      if (usePrec)
        preconditioner = rcp( new MueLu::RefMaxwell<SC,LO,GO,NO>(SM_Matrix,D0_Matrix,M0inv_Matrix,
                                                                 M1_Matrix,nullspace,coords,params) );


#ifdef HAVE_MUELU_TPETRA
      {
        // A test to make sure we can wrap this guy as a MueLu::TpetraOperator
        RCP<Operator> precOp = Teuchos::rcp_dynamic_cast<Operator>(preconditioner);
        MueLu::TpetraOperator<SC,LO,GO,NO> OpT(precOp);
      }
#endif


      // Belos linear problem
#ifdef HAVE_MUELU_BELOS
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(SM_Matrix)); // Turns a Xpetra::Matrix object into a Belos operator

      RCP<Belos::LinearProblem<SC, MV, OP> > problem = rcp( new Belos::LinearProblem<SC, MV, OP>() );
      problem -> setOperator( belosOp );
      Teuchos::RCP<OP> belosPrecOp;
      if (usePrec) {
        belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner)); // Turns a Xpetra::Matrix object into a Belos operator
        problem -> setRightPrec( belosPrecOp );
      }

      problem -> setProblem( X, B );

      bool set = problem->setProblem();
      if (set == false) {
        *out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos solver
      RCP< Belos::SolverManager<SC, MV, OP> > solver;
      RCP< Belos::SolverFactory<SC, MV,OP> > factory = rcp( new  Belos::SolverFactory<SC,MV,OP>() );
      RCP<Teuchos::ParameterList> belosParams
        = rcp( new Teuchos::ParameterList() );
      belosParams->set("Maximum Iterations", 100);
      belosParams->set("Convergence Tolerance",tol);
      belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosParams->set("Output Frequency",1);
      belosParams->set("Output Style",Belos::Brief);
      solver = factory->create(belosSolverType,belosParams);

      comm->barrier();

      auto tm3 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // set problem and solve
      solver -> setProblem( problem );
      Belos::ReturnType status = solver -> solve();
      int iters = solver -> getNumIters();
      success = (iters<50 && status == Belos::Converged);
      if (success)
        *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        *out << "FAILURE! Belos did not converge fast enough." << std::endl;
      tm3 = Teuchos::null;
    }
#ifdef HAVE_MUELU_STRATIMIKOS
    if (solverName == "Stratimikos") {
      auto tm4 = TimeMonitor::getNewTimer("Maxwell: 2 - Build Stratimikos solver");

      // Build the rest of the Stratimikos list
      Teuchos::ParameterList SList;
      SList.set("Linear Solver Type","Belos");
      SList.sublist("Linear Solver Types").sublist("Belos").set("Solver Type", belosSolverType);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Output Frequency",1);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Maximum Iterations",100);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Convergence Tolerance",tol);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Output Style",1);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Verbosity",33);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("VerboseObject").set("Verbosity Level", "medium");
      SList.set("Preconditioner Type","MueLuRefMaxwell");
      params.set("parameterlist: syntax","muelu");
      SList.sublist("Preconditioner Types").set("MueLuRefMaxwell",params);
      // Add matrices to parameterlist
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("D0",D0_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("M0inv",M0inv_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("M1",M1_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("Coordinates",coords);

      // Build Thyra linear algebra objects
      RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(SM_Matrix)->getCrsMatrix());
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(X->getVectorNonConst(0)));
      // TODO: Why do we loose a reference when running this with Epetra?
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(B->getVector(0));

      // Build Stratimikos solver
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
      Stratimikos::enableMueLuRefMaxwell<LocalOrdinal,GlobalOrdinal,Node>(linearSolverBuilder);                // Register MueLu as a Stratimikos preconditioner strategy.
      linearSolverBuilder.setParameterList(rcp(&SList,false));              // Setup solver parameters using a Stratimikos parameter list.

      // Build a new "solver factory" according to the previously specified parameter list.
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

      // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
      Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > thyraInverseA = Thyra::linearOpWithSolve(*solverFactory, thyraA);

      comm->barrier();

      tm4 = Teuchos::null;
      auto tm5 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // Solve Ax = b.
      Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
      std::cout << status << std::endl;

      success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
      tm5 = Teuchos::null;
    }
#endif
    comm->barrier();
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
      if (timingsFormat == "yaml") {
        reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
        reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
      }
      reportParams->set("How to merge timer sets",   "Union");
      reportParams->set("alwaysWriteLocal",          false);
      reportParams->set("writeGlobalStats",          true);
      reportParams->set("writeZeroTimers",           false);
      // FIXME: no "ignoreZeroTimers"

      const std::string filter = "";

      std::ios_base::fmtflags ff(out->flags());
      if (timingsFormat == "table-fixed") *out << std::fixed;
      else * out << std::scientific;
      TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
       *out << std::setiosflags(ff);
    }

    TimeMonitor::clearCounters();
#endif // #ifdef HAVE_MUELU_BELOS
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#else
  return EXIT_SUCCESS;
// #endif // HAVE_TPETRA_INST_INT_INT
#endif
} // main


template<typename Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  return MainWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::main_(clp, lib, argc, argv);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
