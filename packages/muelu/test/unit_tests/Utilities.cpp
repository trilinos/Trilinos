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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_Utilities.hpp>

// This file is intended to house all the tests for MueLu_Utilities.hpp.

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities,MatMatMult_EpetraVsTpetra,Scalar,LocalOrdinal,GlobalOrdinal,Node)
  {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    out << "version: " << MueLu::Version() << std::endl;
    out << "This test compares the matrix matrix multiply between Tpetra and Epetra" << std::endl;

    MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar,GlobalOrdinal,Node);

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    //Calculate result = (Op*Op)*X for Epetra
    GO nx = 37*comm->getSize();
    GO ny = nx;
    RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(nx,ny,Xpetra::UseEpetra);
    RCP<Matrix> OpOp = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Op,false,*Op,false,out);
    RCP<MultiVector> result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    Teuchos::Array<magnitude_type> xnorm(1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Teuchos::Array<magnitude_type> normEpetra(1);
    result->norm2(normEpetra);

    // aid debugging by calculating Op*(Op*X)
    RCP<MultiVector> workVec = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    RCP<MultiVector> check1 = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    Op->apply(*X,*workVec,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Op->apply(*workVec,*check1,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Teuchos::Array<magnitude_type> normCheck1(1);
    check1->norm2(normCheck1);

    //Calculate result = (Op*Op)*X for Tpetra
    Op = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(nx,ny,Xpetra::UseTpetra);
    OpOp = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Op,false,*Op,false,out);
    result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Teuchos::Array<magnitude_type> normTpetra(1);
    result->norm2(normTpetra);

    // aid debugging by calculating Op*(Op*X)
    workVec = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    RCP<MultiVector> check2 = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    Op->apply(*X,*workVec,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Op->apply(*workVec,*check2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    Teuchos::Array<magnitude_type> normCheck2(1);
    check2->norm2(normCheck2);

    TEST_FLOATING_EQUALITY(normEpetra[0], normTpetra[0], 1e-12);
    out << "Epetra ||A*(A*x)|| = " << normCheck1[0] << std::endl;
    out << "Tpetra ||A*(A*x)|| = " << normCheck2[0] << std::endl;
#   else
    out << "Skipping test because some required packages are not enabled (Tpetra, EpetraExt)." << std::endl;
#   endif

  } //EpetraVersusTpetra

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities,DetectDirichletRows,Scalar,LocalOrdinal,GlobalOrdinal,Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(100);
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar>  values;

    LocalOrdinal localRowToZero = 5;
    A->resumeFill();
    A->getLocalRowView(localRowToZero, indices, values);
    Array<Scalar> newvalues(values.size(),TST::zero());
    for (int j = 0; j < indices.size(); j++)
      //keep diagonal
      if (indices[j] == localRowToZero) newvalues[j] = values[j];
    A->replaceLocalValues(localRowToZero,indices,newvalues);

    A->fillComplete();

    ArrayRCP<const bool> drows = Utilities::DetectDirichletRows(*A);
    TEST_EQUALITY(drows[localRowToZero], true);
    TEST_EQUALITY(drows[localRowToZero-1], false);

    A->resumeFill();
    A->getLocalRowView(localRowToZero, indices, values);
    for (int j = 0; j < indices.size(); j++)
      //keep diagonal
      if (indices[j] == localRowToZero) newvalues[j] = values[j];
      else newvalues[j] = Teuchos::as<Scalar>(0.25);
    A->replaceLocalValues(localRowToZero,indices,newvalues);

    //row 5 should not be Dirichlet
    drows = Utilities::DetectDirichletRows(*A,TST::magnitude(0.24));
    TEST_EQUALITY(drows[localRowToZero], false);
    TEST_EQUALITY(drows[localRowToZero-1], false);

    //row 5 should be Dirichlet
    drows = Utilities::DetectDirichletRows(*A,TST::magnitude(0.26));
    TEST_EQUALITY(drows[localRowToZero], true);
    TEST_EQUALITY(drows[localRowToZero-1], false);

  } //DetectDirichletRows

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities,GetDiagonalInverse,Scalar,LocalOrdinal,GlobalOrdinal,Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // blocked diagonal operator (Xpetra)
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);

    RCP<const BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
    TEST_EQUALITY(bA != Teuchos::null, true);

    RCP<Vector> diagInv = Utilities::GetMatrixDiagonalInverse(*bA);
    Teuchos::ArrayRCP<const Scalar> diagInvData = diagInv->getData(0);
    for(size_t i = 0; i < diagInv->getLocalLength(); ++i) {
      if(i >= 0  && i < 5 ) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0),true);
      if(i >= 5  && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5),true);
      if(i >= 10 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0/3.0),true);
    }
    TEST_EQUALITY(diagInv->getMap()->isSameAs(*(bA->getRangeMapExtractor()->getFullMap())),true);

    A = Teuchos::null; bA = Teuchos::null; diagInv = Teuchos::null;

    // blocked diagonal operator (Thyra)
    A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);

    bA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
    TEST_EQUALITY(bA != Teuchos::null, true);

    diagInv = Utilities::GetMatrixDiagonalInverse(*bA);
    diagInvData = diagInv->getData(0);
    for(size_t i = 0; i < diagInv->getLocalLength(); ++i) {
      if(i >= 0  && i < 5 ) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0),true);
      if(i >= 5  && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5),true);
      if(i >= 10 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0/3.0),true);
    }
    TEST_EQUALITY(diagInv->getMap()->isSameAs(*(bA->getRangeMapExtractor()->getFullMap())),true);

    // reordered blocked diagonal operator (Xpetra)
    A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 2 0] 1 ]");
    RCP<const BlockedCrsMatrix> bAA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
    bA = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(buildReorderedBlockedCrsMatrix(brm, bAA));

    TEST_EQUALITY(bA->Rows(),2);
    TEST_EQUALITY(bA->Cols(),2);

    TEST_EQUALITY(bA->getRangeMapExtractor()->getThyraMode(),false);
    TEST_EQUALITY(bA->getDomainMapExtractor()->getThyraMode(),false);

    diagInv = Utilities::GetMatrixDiagonalInverse(*bA);
    diagInvData = diagInv->getData(0);
    for(size_t i = 0; i < diagInv->getLocalLength(); ++i) {
      if(i >= 10 && i < 15) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0),true);
      if(i >= 15 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5),true);
      if(i >= 0  && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0/3.0),true);
    }
    TEST_EQUALITY(diagInv->getMap()->isSameAs(*(bA->getRangeMapExtractor()->getFullMap())),true);


    // reordered blocked diagonal operator (Thyra)
    A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    brm = Xpetra::blockedReorderFromString("[ [ 2 0] 1 ]");
    bAA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
    bA = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(buildReorderedBlockedCrsMatrix(brm, bAA));

    TEST_EQUALITY(bA->Rows(),2);
    TEST_EQUALITY(bA->Cols(),2);

    TEST_EQUALITY(bA->getRangeMapExtractor()->getThyraMode(),true);
    TEST_EQUALITY(bA->getDomainMapExtractor()->getThyraMode(),true);

    diagInv = Utilities::GetMatrixDiagonalInverse(*bA);
    diagInvData = diagInv->getData(0);
    for(size_t i = 0; i < diagInv->getLocalLength(); ++i) {
      if(i >= 10 && i < 15) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0),true);
      if(i >= 15 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5),true);
      if(i >= 0  && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0/3.0),true);
    }
    TEST_EQUALITY(diagInv->getMap()->isSameAs(*(bA->getRangeMapExtractor()->getFullMap())),true);

  } // GetDiagonalInverse

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities,GetLumpedDiagonal,Scalar,LocalOrdinal,GlobalOrdinal,Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    std::vector<RCP<const Map> > maps = std::vector<RCP<const Map> >(3, Teuchos::null);
    maps[0] = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    maps[1] = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    maps[2] = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    RCP<Matrix> A00 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], 4.0, -1.0, -1.0, lib);
    RCP<Matrix> A01 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A10 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A11 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], 4.0, -1.0, -1.0, lib);
    RCP<Matrix> A12 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A21 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A22 = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], 4.0, -1.0, -1.0, lib);

    // create map extractor
    // To generate the Thyra style map extractor we do not need a full map but only the
    // information about the Map details (i.e. lib and indexBase). We can extract this
    // information from maps[0]
    Teuchos::RCP<const MapExtractor > rgMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    Teuchos::RCP<const MapExtractor > doMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    // build blocked operator
    Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor,doMapExtractor,5));
    bop->setMatrix(Teuchos::as<size_t>(0),Teuchos::as<size_t>(0),A00);
    bop->setMatrix(Teuchos::as<size_t>(0),Teuchos::as<size_t>(1),A01);
    bop->setMatrix(Teuchos::as<size_t>(1),Teuchos::as<size_t>(0),A10);
    bop->setMatrix(Teuchos::as<size_t>(1),Teuchos::as<size_t>(1),A11);
    bop->setMatrix(Teuchos::as<size_t>(1),Teuchos::as<size_t>(2),A12);
    bop->setMatrix(Teuchos::as<size_t>(2),Teuchos::as<size_t>(1),A21);
    bop->setMatrix(Teuchos::as<size_t>(2),Teuchos::as<size_t>(2),A22);
    bop->fillComplete();

    // blocked diagonal operator (Thyra)
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(bop);
    TEST_EQUALITY(diagLumped->getMap()->isSameAs(*(bop->getRangeMapExtractor()->getFullMap())),true);

    RCP<Vector> diagLumpedPart0 = bop->getRangeMapExtractor()->ExtractVector(diagLumped,0,false);
    RCP<Vector> diagLumpedPart1 = bop->getRangeMapExtractor()->ExtractVector(diagLumped,1,false);
    RCP<Vector> diagLumpedPart2 = bop->getRangeMapExtractor()->ExtractVector(diagLumped,2,false);
    Teuchos::ArrayRCP<const Scalar> diagLumpedPart0Data = diagLumpedPart0->getData(0);
    Teuchos::ArrayRCP<const Scalar> diagLumpedPart1Data = diagLumpedPart1->getData(0);
    Teuchos::ArrayRCP<const Scalar> diagLumpedPart2Data = diagLumpedPart2->getData(0);
    TEST_EQUALITY(diagLumpedPart0->getLocalLength(),diagLumpedPart1->getLocalLength());
    TEST_EQUALITY(diagLumpedPart0->getLocalLength(),diagLumpedPart2->getLocalLength());
    for(size_t i = 1; i < diagLumpedPart0->getLocalLength() - 1; ++i) {
      TEST_EQUALITY(diagLumpedPart0Data[i],Teuchos::as<Scalar>(7.0));
      TEST_EQUALITY(diagLumpedPart1Data[i],Teuchos::as<Scalar>(8.0));
      TEST_EQUALITY(diagLumpedPart2Data[i],Teuchos::as<Scalar>(7.0));
    }

    LocalOrdinal lastElement = diagLumpedPart0->getLocalLength() - 1;

    if(comm->getSize() == 1) {
      TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(7.0));
      TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(6.0));

      TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(7.0));
      TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(6.0));
    } else {

      if (comm->getRank() == 0) {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(8.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(7.0));

      } else if (comm->getRank() == comm->getSize() - 1) {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(8.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(7.0));

        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(6.0));
      } else {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(8.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(7.0));

        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(8.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(7.0));
      }
    }

    // test reordered operator
    Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 2 0] 1 ]");
    Teuchos::RCP<const BlockedCrsMatrix> bAA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(bop);
    Teuchos::RCP<const BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(buildReorderedBlockedCrsMatrix(brm, bAA));

    TEST_EQUALITY(bA->Rows(),2);
    TEST_EQUALITY(bA->Cols(),2);

    TEST_EQUALITY(bA->getRangeMapExtractor()->getThyraMode(),true);
    TEST_EQUALITY(bA->getDomainMapExtractor()->getThyraMode(),true);

    RCP<Vector> diagLumped2 = Utilities::GetLumpedMatrixDiagonal(bA);


    TEST_EQUALITY(diagLumped2->getMap()->isSameAs(*(bA->getRangeMapExtractor()->getFullMap())),true);

    // the following is only true, since all blocks are the same size and we have
    // Thyra maps?? compare max global ids!
    TEST_EQUALITY(diagLumped2->getMap()->isSameAs(*(diagLumped2->getMap())),true);

    TEST_EQUALITY(diagLumped2->getMap()->getMaxAllGlobalIndex(),comm->getSize() * 300 - 1);
    TEST_EQUALITY(diagLumped2->getMap()->getMinGlobalIndex(),comm->getRank() * 100);

    RCP<const BlockedCrsMatrix> bA0 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(bA->getMatrix(0,0));
    diagLumpedPart0  = bA->getRangeMapExtractor()->ExtractVector(diagLumped2,0,false);
    RCP<Vector> diagLumpedPart00 = bA0->getRangeMapExtractor()->ExtractVector(diagLumpedPart0,0,false);
    RCP<Vector> diagLumpedPart01 = bA0->getRangeMapExtractor()->ExtractVector(diagLumpedPart0,1,false);
    diagLumpedPart1  = bA->getRangeMapExtractor()->ExtractVector(diagLumped2,1,false);

    diagLumpedPart0Data = diagLumpedPart00->getData(0);
    diagLumpedPart1Data = diagLumpedPart01->getData(0);
    diagLumpedPart2Data = diagLumpedPart1->getData(0);
    TEST_EQUALITY(diagLumpedPart0->getLocalLength(),2*diagLumpedPart1->getLocalLength());
    for(size_t i = 1; i < 100 - 1; ++i) {
      TEST_EQUALITY(diagLumpedPart0Data[i],Teuchos::as<Scalar>(7.0));
      TEST_EQUALITY(diagLumpedPart1Data[i],Teuchos::as<Scalar>(7.0));
      TEST_EQUALITY(diagLumpedPart2Data[i],Teuchos::as<Scalar>(8.0));
    }

    lastElement = diagLumpedPart00->getLocalLength() - 1;

    if(comm->getSize() == 1) {
      TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(7.0));

      TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(6.0));
      TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(7.0));
    } else {

      if (comm->getRank() == 0) {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(8.0));

      } else if (comm->getRank() == comm->getSize() - 1) {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(8.0));

        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(6.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(7.0));
      } else {
        TEST_EQUALITY(diagLumpedPart0Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[0],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[0],Teuchos::as<Scalar>(8.0));

        TEST_EQUALITY(diagLumpedPart0Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart1Data[lastElement],Teuchos::as<Scalar>(7.0));
        TEST_EQUALITY(diagLumpedPart2Data[lastElement],Teuchos::as<Scalar>(8.0));
      }
    }
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities,GetInverse,Scalar,LocalOrdinal,GlobalOrdinal,Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    RCP<Map> m  = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(lib, 100, 0, comm);

    RCP<Vector> v  = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(m, true);
    RCP<Vector> tv = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(m, true);
    Teuchos::ArrayRCP<Scalar> vData  = v->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> tvData = tv->getDataNonConst(0);
    for(LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(v->getLocalLength()); ++i) {
      vData[i] = Teuchos::as<Scalar>(i+1);
      tvData[i] = Teuchos::ScalarTraits<Scalar>::one() / Teuchos::as<Scalar>(i+1);
    }

    RCP<Vector> inv = Utilities::GetInverse(v);

    tv->update(1.0,*inv,-1.0);
    TEST_EQUALITY(tv->norm1(),Teuchos::ScalarTraits<Scalar>::zero());
    TEST_EQUALITY(tv->norm2(),Teuchos::ScalarTraits<Scalar>::zero());
    TEST_EQUALITY(tv->normInf(),Teuchos::ScalarTraits<Scalar>::zero());
  }

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
         TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities,MatMatMult_EpetraVsTpetra,Scalar,LO,GO,Node) \
         TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities,DetectDirichletRows,Scalar,LO,GO,Node) \
         TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities,GetDiagonalInverse,Scalar,LO,GO,Node) \
         TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities,GetLumpedDiagonal,Scalar,LO,GO,Node) \
         TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities,GetInverse,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>

}//namespace MueLuTests

