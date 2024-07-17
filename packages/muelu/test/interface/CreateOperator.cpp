// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cstdlib>
#include <fstream>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Tpetra_Operator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <MueLu_CreateEpetraPreconditioner.hpp>
#endif
#include <MueLu_TestHelpers.hpp>

const std::string thickSeparator = "==========================================================================================================================";
const std::string thinSeparator  = "--------------------------------------------------------------------------------------------------------------------------";

const std::string prefSeparator = "=====================================";

namespace MueLuExamples {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void setup_system_list(Xpetra::UnderlyingLib& lib, Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A, Teuchos::ParameterList& mueluList, const std::string& fname) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  int myRank = A->getRowMap()->getComm()->getRank();

  std::filebuf buffer;
  std::streambuf* oldbuffer = NULL;

  typedef Tpetra::Operator<SC, LO, GO, NO> Tpetra_Operator;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> Tpetra_CrsMatrix;
  typedef Tpetra::Vector<SC, LO, GO, NO> Tpetra_Vector;
  typedef Tpetra::MultiVector<SC, LO, GO, NO> Tpetra_MultiVector;
  if (lib == Xpetra::UseTpetra) {
    if (myRank == 0) {
      // Redirect output
      buffer.open((fname + ".out").c_str(), std::ios::out);
      oldbuffer = std::cout.rdbuf(&buffer);
    }

    RCP<Tpetra_CrsMatrix> At = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(A);
    RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > opA(At);
    RCP<Tpetra_Operator> Mt = MueLu::CreateTpetraPreconditioner(opA, mueluList);

    if (myRank == 0) {
      // Redirect output back
      std::cout.rdbuf(oldbuffer);
      buffer.close();
    }
  }
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_SERIAL)
  if (lib == Xpetra::UseEpetra) {
    if (myRank == 0) {
      // Redirect output
      buffer.open((fname + ".out").c_str(), std::ios::out);
      oldbuffer = std::cout.rdbuf(&buffer);
    }

    RCP<Epetra_CrsMatrix> Ae = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(A);
    RCP<Epetra_Operator> Me  = MueLu::CreateEpetraPreconditioner(Ae, mueluList);

    if (myRank == 0) {
      // Redirect output back
      std::cout.rdbuf(oldbuffer);
      buffer.close();
    }
  }
#endif
}

// This routine generate's the user's original A matrix and nullspace
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void generate_user_matrix_and_nullspace(std::string& matrixType, Xpetra::UnderlyingLib& lib, Teuchos::ParameterList& galeriList, Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& nullspace, Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& coordinates) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out       = *fancy;

  RCP<const Map> map;
  if (matrixType == "Laplace1D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian1D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);

  } else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);

  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
  }

  // Expand map to do multiple DOF per node for block problems
  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, (matrixType == "Elasticity2D" ? 2 : 3));

  out << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
      << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
      << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
      << "========================================================" << std::endl;

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, map, galeriList);

  A = Pr->BuildMatrix();

  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
    nullspace = Pr->BuildNullspace();
    A->SetFixedBlockSize((matrixType == "Elasticity2D") ? 2 : 3);
  }
}

void run_sed(const std::string& pattern, const std::string& baseFile) {
  // sed behaviour differs between Mac and Linux
  // You can run "sed -i 's//' " in Linux, but you always have to specify
  // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'
  std::string sed_pref = "sed -i ";
#ifdef __APPLE__
  sed_pref = sed_pref + "\"\" ";
#endif
  int ret_val = 0;
  ret_val     = system((sed_pref + pattern + " " + baseFile + ".gold_filtered").c_str());
  TEUCHOS_ASSERT_EQUALITY(ret_val, 0);
  ret_val = system((sed_pref + pattern + " " + baseFile + ".out_filtered").c_str());
  TEUCHOS_ASSERT_EQUALITY(ret_val, 0);
}

bool compare_to_gold(int myRank, const std::string& baseFile) {
  bool failed = false;
  if (myRank == 0) {
    // Create a copy of outputs
    std::string cmd = "cp -f ";
    int ret_val     = 0;
    ret_val         = system((cmd + baseFile + ".gold " + baseFile + ".gold_filtered").c_str());
    TEUCHOS_ASSERT_EQUALITY(ret_val, 0);
    ret_val = system((cmd + baseFile + ".out " + baseFile + ".out_filtered").c_str());
    TEUCHOS_ASSERT_EQUALITY(ret_val, 0);

    // Tpetra produces different eigenvalues in Chebyshev due to using
    // std::rand() for generating random vectors, which may be initialized
    // using different seed, and may have different algorithm from one
    // gcc version to another, or to anogther compiler (like clang)
    // This leads to us always failing this test.
    // NOTE1 : Epetra, on the other hand, rolls out its out random number
    // generator, which always produces same results

    // make sure complex tests pass
    run_sed("'s/relaxation: damping factor = (1,0)/relaxation: damping factor = 1/'", baseFile);
    run_sed("'s/damping factor: (1,0)/damping factor: 1/'", baseFile);
    run_sed("'s/relaxation: min diagonal value = (0,0)/relaxation: min diagonal value = 0/'", baseFile);

    // Ignore the value of "lambdaMax"
    run_sed("'s/lambdaMax: [0-9]*.[0-9]*/lambdaMax = <ignored>/'", baseFile);

    // Ignore the value of "lambdaMin"
    run_sed("'s/lambdaMin: [0-9]*.[0-9]*/lambdaMin = <ignored>/'", baseFile);

    // Ignore the value of "chebyshev: max eigenvalue"
    // NOTE: we skip lines with default value ([default])
    run_sed("'/[default]/! s/chebyshev: max eigenvalue = [0-9]*.[0-9]*/chebyshev: max eigenvalue = <ignored>/'", baseFile);

    // Ignore the exact type of direct solver (it is selected semi-automatically
    // depending on how Trilinos was configured
    run_sed("'s/Amesos\\([2]*\\)Smoother{type = .*}/Amesos\\1Smoother{type = <ignored>}/'", baseFile);
    run_sed("'s/SuperLU solver interface, direct solve/<Direct> solver interface/'", baseFile);
    run_sed("'s/KLU2 solver interface/<Direct> solver interface/'", baseFile);
    run_sed("'s/Basker solver interface/<Direct> solver interface/'", baseFile);

    // The smoother complexity depends on the coarse solver.
    run_sed("'s/Smoother complexity = [0-9][0-9]*.[0-9]*/Smoother complexity = <ignored>/'", baseFile);

    // Nuke all pointers
    run_sed("'s/0x[0-9a-f]*//g'", baseFile);

    // Strip template args for some classes
    std::vector<std::string> classes;
    classes.push_back("Xpetra::Matrix");
    classes.push_back("MueLu::Constraint");
    classes.push_back("MueLu::SmootherPrototype");
    for (size_t q = 0; q < classes.size(); q++)
      run_sed("'s/" + classes[q] + "<.*>/" + classes[q] + "<ignored> >/'", baseFile);

#ifdef __APPLE__
    // Some Macs print outs ptrs as 0x0 instead of 0, fix that
    run_sed("'/RCP/ s/=0x0/=0/g'", baseFile);
#endif

    // Run comparison (ignoring whitespaces)
    cmd     = "diff -u -w -I\"^\\s*$\" " + baseFile + ".gold_filtered " + baseFile + ".out_filtered";
    int ret = system(cmd.c_str());
    if (ret)
      failed = true;

    std::cout << baseFile << ": " << (ret ? "failed" : "passed") << std::endl;
  }

  return !failed;
}

}  // namespace MueLuExamples

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  bool success = true;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int numProc                         = comm->getSize();
    int myRank                          = comm->getRank();
    RCP<Teuchos::FancyOStream> fancy    = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out          = *fancy;

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
    ::Xpetra::Parameters xpetraParameters(clp);

    bool useKokkos = false;
    if (lib == Xpetra::UseTpetra) {
      useKokkos = !Node::is_serial;
    }
    clp.setOption("useKokkosRefactor", "noKokkosRefactor", &useKokkos, "use kokkos refactor");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    ParameterList galeriList = galeriParameters.GetParameterList();
    out << thickSeparator << std::endl
        << xpetraParameters << galeriParameters;

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<const Map> map;
    RCP<Matrix> A, P, R, Ac;
    RCP<MultiVector> nullspace;
    RCP<RealValuedMultiVector> coordinates0;
    RCP<RealValuedMultiVector> coordinates1;
    std::string matrixType = galeriParameters.GetMatrixType();
    MueLuExamples::generate_user_matrix_and_nullspace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(matrixType, lib, galeriList, comm, A, nullspace, coordinates0);
    map = A->getRowMap();

    std::string prefix;
    if (useKokkos) {
      if (TYPE_EQUAL(Scalar, std::complex<double>) || TYPE_EQUAL(Scalar, std::complex<float>)) {
        prefix = "kokkos-complex/";
      } else {
        prefix = "kokkos/";
      }
    } else {
      if (TYPE_EQUAL(Scalar, std::complex<double>) || TYPE_EQUAL(Scalar, std::complex<float>)) {
        prefix = "complex/";
      } else {
        prefix = "default/";
      }
    }
    std::cout << "Testing folder \"" << prefix << "\"" << std::endl;
    // =========================================================================
    // Solve #1 (standard MueLu)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 1: Standard " << prefSeparator << std::endl;
    {
      std::string fname = prefix + "Output/operator_solve_1_np" + Teuchos::toString(numProc);
      fname             = fname + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

      std::filebuf buffer;
      std::streambuf* oldbuffer = NULL;
      if (myRank == 0) {
        // Redirect output
        buffer.open((fname + ".out").c_str(), std::ios::out);
        oldbuffer = std::cout.rdbuf(&buffer);
      }

      std::srand(12345);

      ParameterList mueluList;
      mueluList.set("verbosity", "interfacetest");
      mueluList.set("coarse: max size", 100);
      mueluList.set("use kokkos refactor", useKokkos);
      mueluList.set("aggregation: deterministic", useKokkos);

      ParameterListInterpreter mueLuFactory(mueluList);
      RCP<Hierarchy> H                              = mueLuFactory.CreateHierarchy();
      Teuchos::RCP<FactoryManagerBase> LevelFactory = mueLuFactory.GetFactoryManager(1);
      H->setlib(lib);
      H->AddNewLevel();
      H->GetLevel(1)->Keep("Nullspace", LevelFactory->GetFactory("Nullspace").get());
      H->GetLevel(1)->Keep("Coordinates", LevelFactory->GetFactory("Coordinates").get());
      H->GetLevel(0)->Set("A", A);
      H->GetLevel(0)->Set("Coordinates", coordinates0);
      mueLuFactory.SetupHierarchy(*H);

      // Extract R, P & Ac for LevelWrap Usage
      H->GetLevel(1)->Get("R", R);
      H->GetLevel(1)->Get("P", P);
      H->GetLevel(1)->Get("A", Ac);
      nullspace    = H->GetLevel(1)->template Get<RCP<MultiVector> >("Nullspace", LevelFactory->GetFactory("Nullspace").get());
      coordinates1 = H->GetLevel(1)->template Get<RCP<RealValuedMultiVector> >("Coordinates", LevelFactory->GetFactory("Coordinates").get());

      if (myRank == 0) {
        // Redirect output back
        std::cout.rdbuf(oldbuffer);
        buffer.close();
      }

      bool passed = MueLuExamples::compare_to_gold(myRank, fname);
      success     = success && passed;
    }

    // =========================================================================
    // Solve #5 (level wrap, the fast way, P, R + Nullspace)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 5: LevelWrap, Fast Way, P, R " << prefSeparator << std::endl;
    {
      std::string fname = prefix + "Output/operator_solve_5_np" + Teuchos::toString(numProc);
      fname             = fname + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

      std::srand(12345);

      ParameterList mueluList;
      mueluList.set("verbosity", "interfacetest");
      mueluList.set("coarse: max size", 100);
      mueluList.set("use kokkos refactor", useKokkos);
      mueluList.set("aggregation: deterministic", useKokkos);
      ParameterList& level0 = mueluList.sublist("level 0");
      level0.set("Coordinates", coordinates0);
      ParameterList& level1 = mueluList.sublist("level 1");
      level1.set("R", R);
      level1.set("P", P);
      level1.set("Nullspace", nullspace);
      level1.set("Coordinates", coordinates1);

      MueLuExamples::setup_system_list<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lib, A, mueluList, fname);

      bool passed = MueLuExamples::compare_to_gold(myRank, fname);
      success     = success && passed;
    }

    // =========================================================================
    // Solve #6 (level wrap, the fast way, P only, explicit transpose)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 6: LevelWrap, Fast Way, P only, explicit transpose " << prefSeparator << std::endl;
    {
      std::string fname = prefix + "Output/operator_solve_6_np" + Teuchos::toString(numProc);
      fname             = fname + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

      std::srand(12345);

      ParameterList mueluList;
      mueluList.set("verbosity", "interfacetest");
      mueluList.set("transpose: use implicit", false);
      mueluList.set("max levels", 4);
      mueluList.set("coarse: max size", 100);
      mueluList.set("use kokkos refactor", useKokkos);
      mueluList.set("aggregation: deterministic", useKokkos);
      ParameterList& level0 = mueluList.sublist("level 0");
      level0.set("Coordinates", coordinates0);
      ParameterList& level1 = mueluList.sublist("level 1");
      level1.set("P", P);
      level1.set("Nullspace", nullspace);
      level1.set("Coordinates", coordinates1);

      MueLuExamples::setup_system_list<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lib, A, mueluList, fname);

      bool passed = MueLuExamples::compare_to_gold(myRank, fname);
      success     = success && passed;
    }

    // =========================================================================
    // Solve #7 (level wrap, the fast way, P only, implicit transpose)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 7: LevelWrap, Fast Way, P only, implicit transpose " << prefSeparator << std::endl;
    {
      std::string fname = prefix + "Output/operator_solve_7_np" + Teuchos::toString(numProc);
      fname             = fname + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

      std::srand(12345);

      ParameterList mueluList;
      mueluList.set("verbosity", "interfacetest");
      mueluList.set("coarse: max size", 100);
      mueluList.set("transpose: use implicit", true);
      mueluList.set("max levels", 2);
      mueluList.set("use kokkos refactor", useKokkos);
      mueluList.set("aggregation: deterministic", useKokkos);
      ParameterList& level0 = mueluList.sublist("level 0");
      level0.set("Coordinates", coordinates0);
      ParameterList& level1 = mueluList.sublist("level 1");
      level1.set("P", P);
      level1.set("Nullspace", nullspace);
      level1.set("Coordinates", coordinates1);

      MueLuExamples::setup_system_list<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lib, A, mueluList, fname);

      bool passed = MueLuExamples::compare_to_gold(myRank, fname);
      success     = success && passed;
    }

    if (myRank == 0)
      std::cout << std::endl
                << "End Result: TEST " << (!success ? "FAILED" : "PASSED") << std::endl;
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
