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

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu_Utilities.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_FactoryManager.hpp>
#include <MueLu_MultiVectorTransferFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_SmootherFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MultiVectorTransferFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Factory> TentativePFact = rcp(new TentativePFactory());
  RCP<Factory> TentativeRFact = rcp(new TransPFactory());  // Use Ptent for coordinate projection

  RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO>("Coordinates"));
  mvtf->SetFactory("R", TentativeRFact);

  TEST_EQUALITY(mvtf != Teuchos::null, true);
}  // Constructor test

//------------------------------------------------------------------------------------------

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MultiVectorTransferFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  out << "Tests the action of the transfer factory on a vector.  In this test, the transfer is the tentative" << std::endl;
  out << "prolongator, and the vector is all ones.  So the norm of the resulting coarse grid vector should be" << std::endl;
  out << "equal to the number of fine degrees of freedom." << std::endl;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;
  typedef typename Xpetra::MultiVector<magnitude_type, LO, GO, NO> RealValuedMultiVector;

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  GO nx         = 199;
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  RCP<MultiVector> fineOnes = MultiVectorFactory::Build(A->getRowMap(), 1);
  fineOnes->putScalar(1.0);
  fineLevel.Set("onesVector", fineOnes);
  fineLevel.Set("Coordinates", coordinates);

  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  RCP<TransPFactory> RFact              = rcp(new TransPFactory());

  RCP<FactoryManager> M = rcp(new FactoryManager());
  M->SetKokkosRefactor(false);
  M->SetFactory("P", TentativePFact);
  M->SetFactory("Ptent", TentativePFact);
  M->SetFactory("R", RFact);
  //    fineLevel.SetFactoryManager(M);
  coarseLevel.SetFactoryManager(M);

  RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO>("onesVector"));
  mvtf->SetFactory("R", RFact);

  coarseLevel.Request("onesVector", mvtf.get());
  coarseLevel.Request("R", RFact.get());
  coarseLevel.Request("P", TentativePFact.get());

  mvtf->Build(fineLevel, coarseLevel);

  RCP<MultiVector> coarseOnes = coarseLevel.Get<RCP<MultiVector> >("onesVector", mvtf.get());
  Teuchos::Array<magnitude_type> vn(1);
  coarseOnes->norm2(vn);
  TEST_FLOATING_EQUALITY(vn[0] * vn[0], (Teuchos::as<magnitude_type>(fineOnes->getGlobalLength())), 1e-12);
}  // Build test

//------------------------------------------------------------------------------------------

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MultiVectorTransferFactory, ThreeLevels, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
#if !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Ifpack");
#endif
#if !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Ifpack2");
#endif
  out << "version: " << MueLu::Version() << std::endl;

  out << "Tests usage on a three-level hierarchy." << std::endl;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  GO nx         = 199;
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Set up three level hierarchy.
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  RCP<Level> fineLevel = H->GetLevel();
  fineLevel->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  fineLevel->Set("A", A);  // set fine level matrix
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullSpace->putScalar((SC)1.0);
  fineLevel->Set("Nullspace", nullSpace);      // set null space information for finest level
  fineLevel->Set("Coordinates", coordinates);  // set coordinates on finest level

  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  UncoupledAggFact->SetMinNodesPerAggregate(3);
  UncoupledAggFact->SetMaxNeighAlreadySelected(0);
  UncoupledAggFact->SetOrdering("natural");

  RCP<TentativePFactory> PFact = rcp(new TentativePFactory());  // just using plain aggregation
  RCP<Factory> RFact           = rcp(new TransPFactory());
  RCP<RAPFactory> AcFact       = rcp(new RAPFactory());
  H->SetMaxCoarseSize(1);

  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LO)1);
  smootherParamList.set("relaxation: damping factor", (SC)1.0);
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = Teuchos::null;
  AcFact->setVerbLevel(Teuchos::VERB_HIGH);

  FactoryManager M;
  M.SetKokkosRefactor(false);
  M.SetFactory("Aggregates", UncoupledAggFact);
  M.SetFactory("P", PFact);
  M.SetFactory("Ptent", PFact);  // for nullspace
  M.SetFactory("R", RFact);
  M.SetFactory("A", AcFact);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", SmooFact);  // This line avoid dependency to Amesos/Amesos2 for this test.

  // set up the transfer factory
  RCP<MultiVector> fineOnes = MultiVectorFactory::Build(A->getRowMap(), 1);
  fineOnes->putScalar(1.0);
  fineLevel->Set("onesVector", fineOnes);
  RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO>("onesVector"));
  mvtf->SetFactory("R", RFact);
  M.SetFactory("onesVector", mvtf);
  AcFact->AddTransferFactory(mvtf);

  int maxLevels = 3;
  H->Setup(M, 0, maxLevels);

  /*
  //FIXME we probably need to do some requests....
  coarseLevel.Request("onesVector",mvtf.get());
  coarseLevel.Request("R",RFact.get());
  coarseLevel.Request("P",TentativePFact.get());
  */

  /*
     RCP<MultiVector> coarseOnes = coarseLevel.Get<RCP<MultiVector> >("onesVector",mvtf.get());
     Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> vn(1);
     coarseOnes->norm2(vn);

     TEST_FLOATING_EQUALITY(vn[0]*vn[0],((SC)fineOnes->getGlobalLength()),1e-12);
     */
}  // ThreeLevels

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MultiVectorTransferFactory, Constructor, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MultiVectorTransferFactory, Build, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MultiVectorTransferFactory, ThreeLevels, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
