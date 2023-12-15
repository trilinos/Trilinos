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
#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_NullspaceFactory.hpp"
//#include "MueLu_FactoryManager.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"
#include "MueLu_LineDetectionFactory.hpp"
#include "MueLu_SemiCoarsenPFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SemiCoarsenPFactoryWithSemiRestriction, TestSemiRestrict, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  // Don't test for complex - matrix reader won't work
  if (TST::isComplex) {
    success = true;
    return;
  }

  // MAKE SURE lines are not split across processors !!
  //    ==> nMatrixRows/(nProcs*lineLength*blkSize) must be an integer
  //    ==> use vertical orientation
  //    GO lineLength = 7, blkSize = 1, nMatrixRows= 84;
  int lineLength = 7, blkSize = 2, nMatrixRows = 168;

  // Set up factories for test

  RCP<LineDetectionFactory> LineDetectionFact = rcp(new LineDetectionFactory());
  LineDetectionFact->SetParameter("linedetection: orientation",
                                  Teuchos::ParameterEntry(std::string("vertical")));
  LineDetectionFact->SetParameter("linedetection: num layers",
                                  Teuchos::ParameterEntry(lineLength));
  LineDetectionFact->SetVerbLevel(MueLu::Extreme);
  RCP<SemiCoarsenPFactory> SemiCoarsenPFact = rcp(new SemiCoarsenPFactory());
  SemiCoarsenPFact->SetVerbLevel(MueLu::Extreme);
  SemiCoarsenPFact->SetParameter("semicoarsen: coarsen rate", Teuchos::ParameterEntry(2));
  SemiCoarsenPFact->SetParameter("semicoarsen: calculate nonsym restriction", Teuchos::ParameterEntry(true));
  SemiCoarsenPFact->SetFactory("LineDetection_VertLineIds", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("LineDetection_Layers", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("CoarseNumZLayers", LineDetectionFact);
  // seem to need two factories to avoid failing multipleCallCheck on some machines?
  RCP<LineDetectionFactory> LineDetectionFact2 = rcp(new LineDetectionFactory());
  LineDetectionFact2->SetParameter("linedetection: orientation",
                                   Teuchos::ParameterEntry(std::string("vertical")));
  LineDetectionFact2->SetParameter("linedetection: num layers",
                                   Teuchos::ParameterEntry(lineLength));
  LineDetectionFact2->SetVerbLevel(MueLu::Extreme);
  RCP<SemiCoarsenPFactory> SemiCoarsenPFact2 = rcp(new SemiCoarsenPFactory());
  SemiCoarsenPFact2->SetVerbLevel(MueLu::Extreme);
  SemiCoarsenPFact2->SetParameter("semicoarsen: coarsen rate", Teuchos::ParameterEntry(2));
  SemiCoarsenPFact2->SetParameter("semicoarsen: calculate nonsym restriction", Teuchos::ParameterEntry(true));
  SemiCoarsenPFact2->SetFactory("LineDetection_VertLineIds", LineDetectionFact2);
  SemiCoarsenPFact2->SetFactory("LineDetection_Layers", LineDetectionFact2);
  SemiCoarsenPFact2->SetFactory("CoarseNumZLayers", LineDetectionFact2);

  const RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), nMatrixRows, 0, comm);
  //    RCP<Matrix> Op =  Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/semiRtest2.mm", map);
  RCP<Matrix> Op = Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/semiRblkTestver.mm", map);
  Op->SetFixedBlockSize(blkSize);
  RCP<Matrix> transposedOp = Utilities::Transpose(*Op);

  SC one                     = Teuchos::ScalarTraits<SC>::one();
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(Op->getRowMap(), blkSize);

  for (int i = 0; i < blkSize; i++) {
    RCP<const Map> domainMap = Op->getDomainMap();
    GO indexBase             = domainMap->getIndexBase();

    ArrayRCP<SC> nsData = nullSpace->getDataNonConst(i);
    for (int j = 0; j < nsData.size(); j++) {
      GO GID = domainMap->getGlobalElement(j) - indexBase;
      if ((GID - i) % blkSize == 0) nsData[j] = one;
    }
  }

  Level fineLevel1, coarseLevel1, fineLevel2, coarseLevel2;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel1, coarseLevel1);
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel2, coarseLevel2);
  fineLevel1.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel1.SetFactoryManager(Teuchos::null);
  fineLevel1.Set("A", Op);
  fineLevel1.Set("Nullspace", nullSpace);
  coarseLevel1.Request("P", SemiCoarsenPFact.get());
  SemiCoarsenPFact->Build(fineLevel1, coarseLevel1);

  fineLevel2.Set("A", transposedOp);
  fineLevel2.Set("Nullspace", nullSpace);
  fineLevel2.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel2.SetFactoryManager(Teuchos::null);
  coarseLevel2.Request("P", SemiCoarsenPFact2.get());
  coarseLevel2.Request("RfromPfactory", SemiCoarsenPFact2.get());
  SemiCoarsenPFact2->Build(fineLevel2, coarseLevel2);

  RCP<Matrix> firstP            = coarseLevel1.Get<RCP<Matrix> >("P", &(*SemiCoarsenPFact));
  RCP<Matrix> secondR           = coarseLevel2.Get<RCP<Matrix> >("RfromPfactory", &(*SemiCoarsenPFact2));
  RCP<Matrix> transposedSecondR = Utilities::Transpose(*secondR);
  // Xpetra::IO<SC,LO,GO,NO>::Write("firstP.mm",*firstP);
  // Xpetra::IO<SC,LO,GO,NO>::Write("secondRtrans.mm",*transposedSecondR);

  Utilities::SetRandomSeed(*comm);
  RCP<MultiVector> workVec1 = MultiVectorFactory::Build(firstP->getRangeMap(), 1);
  RCP<MultiVector> workVec2 = MultiVectorFactory::Build(transposedSecondR->getRangeMap(), 1);
  RCP<MultiVector> X        = MultiVectorFactory::Build(firstP->getDomainMap(), 1);
  X->randomize();
  firstP->apply(*X, *workVec1, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  transposedSecondR->apply(*X, *workVec2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  RCP<MultiVector> diff = MultiVectorFactory::Build(firstP->getRowMap(), 1);
  diff->putScalar(0.0);
  // diff = workVec1 + (-1.0)*(workVec2) + 0*diff
  diff->update(1.0, *workVec1, -1.0, *workVec2, 0.0);
  Teuchos::Array<typename TST::magnitudeType> norms(1);
  diff->norm2(norms);
  out << "||diff"
      << "||_2 = " << norms[0] << std::endl;
  TEST_EQUALITY(norms[0] < 1000 * TMT::eps(), true);

}  // NullSpace test

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SemiCoarsenPFactoryWithSemiRestriction, TestSemiRestrict, Scalar, LO, GO, Node)
#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
