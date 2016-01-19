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

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

// This example demonstrates how to reuse some parts of a classical SA multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.

#ifdef HAVE_MUELU_TPETRA
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_CreateEpetraPreconditioner.hpp>
#endif

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
CreateHierarchy(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A, Teuchos::ParameterList& paramList, Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > coords = Teuchos::null) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;

  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

  if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >     tA = Utilities::Op2NonConstTpetraCrs(A);
    RCP<MueLu::TpetraOperator<SC, LO, GO, NO> > tH = MueLu::CreateTpetraPreconditioner(tA, paramList, Utilities::MV2NonConstTpetraMV(coords));

    return tH->GetHierarchy();
#else
    throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif // HAVE_MUELU_TPETRA
  }

  XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
  XPETRA_FACTORY_END;
}

#ifdef HAVE_MUELU_EPETRA
template<>
Teuchos::RCP<MueLu::Hierarchy<double,int,int,Xpetra::EpetraNode> >
CreateHierarchy(Teuchos::RCP<Xpetra::Matrix<double,int,int,Xpetra::EpetraNode> > A, Teuchos::ParameterList& paramList, Teuchos::RCP<Xpetra::MultiVector<double,int,int,Xpetra::EpetraNode> > coords) {
  typedef double                Scalar;
  typedef int                   LocalOrdinal;
  typedef int                   GlobalOrdinal;
  typedef Xpetra::EpetraNode    Node;
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;

  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

  if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >     tA = Utilities::Op2NonConstTpetraCrs(A);
    RCP<MueLu::TpetraOperator<SC, LO, GO, NO> > tH = MueLu::CreateTpetraPreconditioner(tA, paramList, Utilities::MV2NonConstTpetraMV(coords));

    return tH->GetHierarchy();
#else
    throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
  }

  if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
    RCP<Epetra_CrsMatrix>      eA = Utilities::Op2NonConstEpetraCrs(A);
    RCP<MueLu::EpetraOperator> eH = MueLu::CreateEpetraPreconditioner(eA, paramList, Utilities::MV2NonConstEpetraMV(coords));
    return eH->GetHierarchy();
#else
    throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
  }

  XPETRA_FACTORY_END;
}
#endif // HAVE_MUELU_EPETRA

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
ReuseHierarchy(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A, MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>& H) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;

  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

  if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
    MueLu::TpetraOperator<SC, LO, GO, NO> tH(Teuchos::rcpFromRef(H));
    MueLu::ReuseTpetraPreconditioner(Utilities::Op2NonConstTpetraCrs(A), tH);
    return;
#else
    throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
  }

  XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
  XPETRA_FACTORY_END;
}

#ifdef HAVE_MUELU_EPETRA
template<>
void
ReuseHierarchy(Teuchos::RCP<Xpetra::Matrix<double,int,int,Xpetra::EpetraNode> > A, MueLu::Hierarchy<double,int,int,Xpetra::EpetraNode>& H) {
  typedef double                Scalar;
  typedef int                   LocalOrdinal;
  typedef int                   GlobalOrdinal;
  typedef Xpetra::EpetraNode    Node;
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;

  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

  if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
    MueLu::TpetraOperator<SC, LO, GO, NO> tH(Teuchos::rcpFromRef(H));
    MueLu::ReuseTpetraPreconditioner(Utilities::Op2NonConstTpetraCrs(A), tH);
    return;
#else
    throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
  }

  if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
    MueLu::EpetraOperator eH(Teuchos::rcpFromRef(H));
    MueLu::ReuseEpetraPreconditioner(Utilities::Op2NonConstEpetraCrs(A), eH);
    return;
#else
    throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
  }

  XPETRA_FACTORY_END;
}
#endif

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC one = STS::one(), zero = STS::zero();

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = ""; clp.setOption("xml",        &xmlFileName, "read parameters from a file");
  int         numRebuilds = 0;  clp.setOption("rebuild",    &numRebuilds, "#times to rebuild hierarchy");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  ParameterList paramList;
  if (xmlFileName != "")
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  // For comments, see Driver.cpp
  out << "========================================================\n" << xpetraParameters << galeriParameters;

  RCP<Matrix>      A;
  RCP<const Map>   map;
  RCP<MultiVector> coordinates;
  RCP<MultiVector> nullspace;

  {
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

    out << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
        << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
        << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
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
  }

  // =========================================================================
  // Setups and solves
  // =========================================================================
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);
  B->setSeed(846930886);
  B->randomize();

  const int nIts = 9;

  std::string thickSeparator = "==========================================================================================================================";
  std::string thinSeparator  = "--------------------------------------------------------------------------------------------------------------------------";

  RCP<TimeMonitor> tm;

  // =========================================================================
  // Solve #1 (no reuse)
  // =========================================================================
  out << thickSeparator << std::endl;
  {
    RCP<Hierarchy> H;

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #1: no reuse")));
    for (int i = 0; i <= numRebuilds; i++) {
      A->SetMaxEigenvalueEstimate(-one);
      H = CreateHierarchy(A, paramList, coordinates);
    }
    tm = Teuchos::null;

    X->putScalar(zero);
    H->Iterate(*B, *X, nIts);
    out << "||Residual|| = " << Utilities::ResidualNorm(*A, *X, *B)[0] << std::endl;
  }

  // =========================================================================
  // Solve #2 (reuse tentative P)
  // =========================================================================
  out << thickSeparator << std::endl;
  {
    A->SetMaxEigenvalueEstimate(-one);

    paramList.set("reuse: type", "tP");

    RCP<Hierarchy> H = CreateHierarchy(A, paramList, coordinates);

    out << thinSeparator << "\nPreconditioner status:" << std::endl;
    H->print(out, MueLu::Extreme);

    // Reuse setup
    RCP<Matrix> Acopy = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(A);

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #2: reuse tentative P")));
    for (int i = 0; i <= numRebuilds; i++) {
      Acopy->SetMaxEigenvalueEstimate(-one);
      ReuseHierarchy(Acopy, *H);
    }
    tm = Teuchos::null;

    X->putScalar(zero);
    H->Iterate(*B, *X, nIts);
    out << "||Residual|| = " << Utilities::ResidualNorm(*A, *X, *B)[0] << std::endl;
  }

  // =========================================================================
  // Solve #3 (reuse smoothed P and R)
  // =========================================================================
  out << thickSeparator << std::endl;
  {
    A->SetMaxEigenvalueEstimate(-one);

    paramList.set("reuse: type", "RP");

    RCP<Hierarchy> H = CreateHierarchy(A, paramList, coordinates);

    out << thinSeparator << "\nPreconditioner status:" << std::endl;
    H->print(out, MueLu::Extreme);

    // Reuse setup
    RCP<Matrix> Acopy = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(A);
    Acopy->SetMaxEigenvalueEstimate(-one);

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #3: reuse smoothed P and R")));
    for (int i = 0; i <= numRebuilds; i++) {
      Acopy->SetMaxEigenvalueEstimate(-one);
      ReuseHierarchy(Acopy, *H);
    }
    tm = Teuchos::null;

    X->putScalar(zero);
    H->Iterate(*B, *X, nIts);
    out << "||Residual|| = " << Utilities::ResidualNorm(*A, *X, *B)[0] << std::endl;
  }

  {
    const bool alwaysWriteLocal = false;
    const bool writeGlobalStats = true;
    const bool writeZeroTimers  = false;
    const bool ignoreZeroTimers = true;
    const std::string filter    = "Solve #";
    TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                           writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
  }

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  bool verbose = true;
  int return_code = EXIT_FAILURE;

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return_code = main_<double,int,int,Xpetra::EpetraNode>(clp, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
      typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
      return_code = main_<double,int,long,Node>(clp, argc, argv);
#else
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
      return_code = main_<double,int,int,Node> (clp, argc, argv);
#elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
      return_code = main_<double,int,long,Node>(clp, argc, argv);
#elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
      return_code = main_<double,int,long long,Node>(clp, argc, argv);
#else
      throw std::runtime_error("Found no suitable instantiation");
#endif
#endif
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, return_code);

  return return_code;
}
