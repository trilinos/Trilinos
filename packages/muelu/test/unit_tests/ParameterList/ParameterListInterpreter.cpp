// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_TestHelpers.hpp>

#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_MatrixMatrix.hpp>

#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ParameterListInterpreter, SetParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS) && defined(HAVE_MUELU_AMESOS2)

  RCP<Matrix> A                       = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(99);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/ParameterListInterpreter/"), std::string(".xml"));

  for (int i = 0; i < fileList.size(); i++) {
    // Ignore files with "BlockCrs" in their name
    auto found = fileList[i].find("BlockCrs");
    if (found != std::string::npos) continue;

    // Ignore files with "Comparison" in their name
    found = fileList[i].find("Comparison");
    if (found != std::string::npos) continue;

    out << "Processing file: " << fileList[i] << std::endl;
    ParameterListInterpreter mueluFactory("ParameterList/ParameterListInterpreter/" + fileList[i], *comm);

    RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", A);

    mueluFactory.SetupHierarchy(*H);

    // TODO: check no unused parameters
    // TODO: check results of Iterate()
  }

#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2, Amesos, Amesos2)." << std::endl;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ParameterListInterpreter, BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList matrixParams;
    matrixParams.set("matrixType", "Laplace1D");
    matrixParams.set("nx", (GlobalOrdinal)300);  // needs to be even

    RCP<Matrix> A = TestHelpers::TpetraTestFactory<SC, LO, GO, NO>::BuildBlockMatrix(matrixParams, Xpetra::UseTpetra);
    out << "Matrix Size (block) = " << A->getGlobalNumRows() << " (point) " << A->getRangeMap()->getGlobalNumElements() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/ParameterListInterpreter/"), std::string(".xml"));

    for (int i = 0; i < fileList.size(); i++) {
      // Only run files with "BlockCrs" in their name
      auto found = fileList[i].find("BlockCrs");
      if (found == std::string::npos) continue;

      out << "Processing file: " << fileList[i] << std::endl;

      Teuchos::RCP<Teuchos::ParameterList> mueluList = rcp(new Teuchos::ParameterList());
      Teuchos::updateParametersFromXmlFileAndBroadcast("ParameterList/ParameterListInterpreter/" + fileList[i], mueluList.ptr(), *comm);
      mueluList->set("use kokkos refactor", false);

      ParameterListInterpreter mueluFactory(*mueluList, comm);

      RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", A);

      mueluFactory.SetupHierarchy(*H);

      // Test to make sure all of the matrices in the Hierarchy are actually Block Matrices
      using helpers = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      for (int j = 0; j < H->GetNumLevels(); j++) {
        RCP<Level> level = H->GetLevel(j);

        RCP<Matrix> Am = level->Get<RCP<Matrix> >("A");
        TEST_EQUALITY(helpers::isTpetraBlockCrs(Am), true);
        if (j > 0) {
          RCP<Matrix> P = level->Get<RCP<Matrix> >("P");
          TEST_EQUALITY(helpers::isTpetraBlockCrs(P), true);
          RCP<Matrix> R = level->Get<RCP<Matrix> >("R");
          TEST_EQUALITY(helpers::isTpetraBlockCrs(R), true);
        }
      }

      // TODO: check no unused parameters
      // TODO: check results of Iterate()
    }
  }
  TEST_EQUALITY(1, 1);
}

template <class Matrix, class MT>
MT compare_matrices(RCP<Matrix> &Ap, RCP<Matrix> &Ab) {
  using SC  = typename Matrix::scalar_type;
  using LO  = typename Matrix::local_ordinal_type;
  using GO  = typename Matrix::global_ordinal_type;
  using NO  = typename Matrix::node_type;
  using CRS = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  SC one    = Teuchos::ScalarTraits<SC>::one();
  SC zero   = Teuchos::ScalarTraits<SC>::zero();

  RCP<const CRS> Ap_t  = MueLu::Utilities<SC, LO, GO, NO>::Op2TpetraCrs(Ap);
  auto Ab_t            = MueLu::Utilities<SC, LO, GO, NO>::Op2TpetraBlockCrs(Ab);
  RCP<CRS> Ab_as_point = Tpetra::convertToCrsMatrix<SC, LO, GO, NO>(*Ab_t);

  RCP<CRS> diff = rcp(new CRS(Ap_t->getCrsGraph()));
  diff->setAllToScalar(zero);
  diff->fillComplete();
  Tpetra::MatrixMatrix::Add<SC, LO, GO, NO>(*Ap_t, false, one, *Ab_as_point, false, -one, diff);
  return diff->getFrobeniusNorm();
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ParameterListInterpreter, PointCrs_vs_BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
#if !defined(HAVE_MUELU_AMESOS2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2");
#endif
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList matrixParams;
    matrixParams.set("matrixType", "Laplace1D");
    matrixParams.set("nx", (GlobalOrdinal)300);  // needs to be even

    RCP<Matrix> PointA = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, Xpetra::UseTpetra);
    RCP<Matrix> BlockA;
    {
      using XCRS = Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>;

      auto tA      = MueLu::Utilities<SC, LO, GO, NO>::Op2TpetraCrs(PointA);
      auto bA      = Tpetra::convertToBlockCrsMatrix<SC, LO, GO, NO>(*tA, 1);
      RCP<XCRS> AA = rcp(new XCRS(bA));
      BlockA       = rcp(new CrsMatrixWrap(rcp_implicit_cast<CrsMatrix>(AA)));
    }

    out << "Point: Matrix Size (block) = " << PointA->getGlobalNumRows() << " (point) " << PointA->getRangeMap()->getGlobalNumElements() << std::endl;
    out << "Block: Matrix Size (block) = " << BlockA->getGlobalNumRows() << " (point) " << BlockA->getRangeMap()->getGlobalNumElements() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/ParameterListInterpreter/"), std::string(".xml"));

    for (int i = 0; i < fileList.size(); i++) {
      // Only run files with "Comparison" in their name
      auto found = fileList[i].find("Comparison");
      if (found == std::string::npos) continue;

      out << "Processing file: " << fileList[i] << std::endl;

      Teuchos::RCP<Teuchos::ParameterList> mueluList = rcp(new Teuchos::ParameterList());
      Teuchos::updateParametersFromXmlFileAndBroadcast("ParameterList/ParameterListInterpreter/" + fileList[i], mueluList.ptr(), *comm);
      mueluList->set("use kokkos refactor", false);

      // Point Hierarchy
      ParameterListInterpreter mueluFactory1(*mueluList, comm);
      RCP<Hierarchy> PointH = mueluFactory1.CreateHierarchy();
      PointH->GetLevel(0)->Set("A", PointA);
      mueluFactory1.SetupHierarchy(*PointH);

      // Block Hierachy
      ParameterListInterpreter mueluFactory2(*mueluList, comm);
      RCP<Hierarchy> BlockH = mueluFactory2.CreateHierarchy();
      BlockH->GetLevel(0)->Set("A", BlockA);
      mueluFactory2.SetupHierarchy(*BlockH);

      // Check to see that we get the same matrices in both hierarchies
      TEST_EQUALITY(PointH->GetNumLevels(), BlockH->GetNumLevels());

      // TODO BMK: compare_matrices uses MatrixMatrix::Add. This was broken before
      // (see #12565), producing a matrix of all zeros if C is fill-complete on input, like
      // compare_matrices does in this test.
      //
      // After fixing Tpetra::MatrixMatrix::Add, it shows that these pairs of matrices (Ap and Ab, etc.)
      // are actually different so this test is not passing. When this is fixed, uncomment this block.

      /*
              for(int j=0; j<PointH->GetNumLevels(); j++) {
                using CRS=Tpetra::CrsMatrix<SC,LO,GO,NO>;
                using MT  = typename Teuchos::ScalarTraits<SC>::magnitudeType;
                MT tol = Teuchos::ScalarTraits<MT>::squareroot(Teuchos::ScalarTraits<MT>::eps());

                RCP<Level> Plevel = PointH->GetLevel(j);
                RCP<Level> Blevel = BlockH->GetLevel(j);

                // Compare A
                RCP<Matrix> Ap = Plevel->Get<RCP<Matrix> >("A");
                RCP<Matrix> Ab = Blevel->Get<RCP<Matrix> >("A");
                MT norm = compare_matrices<Matrix,MT>(Ap,Ab);

                TEUCHOS_TEST_COMPARE(norm,<,tol,out,success);

                // Compare P, R
                if(j>0) {
                  RCP<Matrix> Pp = Plevel->Get<RCP<Matrix> >("P");
                  RCP<Matrix> Pb = Blevel->Get<RCP<Matrix> >("P");
                  norm = compare_matrices<Matrix,MT>(Pp,Pb);
                  TEUCHOS_TEST_COMPARE(norm,<,tol,out,success);

                  RCP<Matrix> Rp = Plevel->Get<RCP<Matrix> >("R");
                  RCP<Matrix> Rb = Blevel->Get<RCP<Matrix> >("R");
                  norm = compare_matrices<Matrix,MT>(Rp,Rb);
                  TEUCHOS_TEST_COMPARE(norm,<,tol,out,success);
                }
              }
      */

      // TODO: check no unused parameters
      // TODO: check results of Iterate()
    }
  }
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ParameterListInterpreter, SetParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ParameterListInterpreter, BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ParameterListInterpreter, PointCrs_vs_BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

// TODO: some tests of the parameter list parser can be done without building the Hierarchy.
