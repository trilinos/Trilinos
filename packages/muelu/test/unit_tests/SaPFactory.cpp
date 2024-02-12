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
#include <Teuchos_DefaultComm.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_SaPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, Test0, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  TEST_EQUALITY(sapFactory != Teuchos::null, true);

  out << *sapFactory << std::endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, EpetraVsTpetra, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2)
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  out << "version: " << MueLu::Version() << std::endl;
  out << "Compare results of Epetra and Tpetra" << std::endl;
  out << "for 3 level AMG solver using smoothed aggregation with" << std::endl;
  out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar, GlobalOrdinal, Node);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::Array<magnitude_type> results(2);

  // run test only on 1 proc
  if (comm->getSize() == 1) {
    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // run Epetra and Tpetra test
    for (int run = 0; run < 2; run++)  // TODO: create a subfunction instead or Tuple of UnderlyingLib
    {
      if (run == 0)
        lib = Xpetra::UseEpetra;
      else
        lib = Xpetra::UseTpetra;

      // generate problem
      LocalOrdinal maxLevels   = 3;
      LocalOrdinal its         = 10;
      GlobalOrdinal nEle       = 63;
      const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
      Teuchos::ParameterList matrixParameters;
      matrixParameters.set("nx", nEle);

      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace1D", map, matrixParameters);
      RCP<Matrix> Op = Pr->BuildMatrix();

      // build nullspace
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
      nullSpace->putScalar((Scalar)1.0);
      Teuchos::Array<magnitude_type> norms(1);
      nullSpace->norm1(norms);
      if (comm->getRank() == 0)
        out << "||NS|| = " << norms[0] << std::endl;

      // fill hierarchy
      RCP<Hierarchy> H = rcp(new Hierarchy());
      H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

      RCP<Level> Finest = H->GetLevel();
      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set("A", Op);                 // set fine level matrix
      Finest->Set("Nullspace", nullSpace);  // set null space information for finest level

      // define transfer operators
      RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
      UncoupledAggFact->SetMinNodesPerAggregate(3);
      UncoupledAggFact->SetMaxNeighAlreadySelected(0);
      UncoupledAggFact->SetOrdering("natural");

      RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
      RCP<SaPFactory> Pfact            = rcp(new SaPFactory());
      RCP<Factory> Rfact               = rcp(new TransPFactory());
      RCP<RAPFactory> Acfact           = rcp(new RAPFactory());
      H->SetMaxCoarseSize(1);

      // setup smoothers
      Teuchos::ParameterList smootherParamList;
      smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
      smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
      smootherParamList.set("relaxation: damping factor", (Scalar)1.0);
      RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
      RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
      Acfact->setVerbLevel(Teuchos::VERB_HIGH);

      RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

      FactoryManager M;
      M.SetKokkosRefactor(false);
      M.SetFactory("P", Pfact);
      M.SetFactory("R", Rfact);
      M.SetFactory("A", Acfact);
      M.SetFactory("Ptent", Ptentfact);
      M.SetFactory("Aggregates", UncoupledAggFact);
      M.SetFactory("Smoother", SmooFact);
      M.SetFactory("CoarseSolver", coarseSolveFact);

      H->Setup(M, 0, maxLevels);

      // test some basic multigrid data
      RCP<Level> coarseLevel = H->GetLevel(1);
      TEST_EQUALITY(coarseLevel->IsRequested("A", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("P", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("R", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("A", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel->IsAvailable("P", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel->IsAvailable("R", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("A", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("P", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("R", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel->IsRequested("P", Pfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("P", Ptentfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("R", Rfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsRequested("A", Acfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("P", Pfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("P", Ptentfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("R", Rfact.get()), false);
      TEST_EQUALITY(coarseLevel->IsAvailable("A", Acfact.get()), false);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("P", Pfact.get()), 0);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("P", Ptentfact.get()), 0);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", SmooFact.get()), 0);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother", SmooFact.get()), 0);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("R", Rfact.get()), 0);
      TEST_EQUALITY(coarseLevel->GetKeepFlag("A", Acfact.get()), 0);
      RCP<Matrix> P1 = coarseLevel->Get<RCP<Matrix> >("P");
      RCP<Matrix> R1 = coarseLevel->Get<RCP<Matrix> >("R");
      TEST_EQUALITY(P1->getGlobalNumRows(), 63);
      TEST_EQUALITY(P1->getGlobalNumCols(), 21);
      TEST_EQUALITY(R1->getGlobalNumRows(), 21);
      TEST_EQUALITY(R1->getGlobalNumCols(), 63);
      RCP<Level> coarseLevel2 = H->GetLevel(2);
      TEST_EQUALITY(coarseLevel2->IsRequested("A", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("P", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("R", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("A", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel2->IsAvailable("P", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("R", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("A", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("P", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), 0);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("R", MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(coarseLevel2->IsRequested("P", Pfact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("P", Ptentfact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsRequested("R", Rfact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("P", Pfact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("P", Ptentfact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother", SmooFact.get()), false);
      TEST_EQUALITY(coarseLevel2->IsAvailable("R", Rfact.get()), false);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("P", Pfact.get()), 0);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("P", Ptentfact.get()), 0);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother", SmooFact.get()), 0);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother", SmooFact.get()), 0);
      TEST_EQUALITY(coarseLevel2->GetKeepFlag("R", Rfact.get()), 0);
      RCP<Matrix> P2 = coarseLevel2->Get<RCP<Matrix> >("P");
      RCP<Matrix> R2 = coarseLevel2->Get<RCP<Matrix> >("R");
      TEST_EQUALITY(P2->getGlobalNumRows(), 21);
      TEST_EQUALITY(P2->getGlobalNumCols(), 7);
      TEST_EQUALITY(R2->getGlobalNumRows(), 7);
      TEST_EQUALITY(R2->getGlobalNumCols(), 21);

      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P1, true, *P1, false, out);

      if (PtentTPtent->haveGlobalConstants()) TEST_EQUALITY(PtentTPtent->getGlobalMaxNumRowEntries() - 3 < 1e-12, true);
      if (P1->haveGlobalConstants()) TEST_EQUALITY(P1->getGlobalMaxNumRowEntries() - 2 < 1e-12, true);
      if (P2->haveGlobalConstants()) TEST_EQUALITY(P2->getGlobalMaxNumRowEntries() - 2 < 1e-12, true);

      // Define RHS
      RCP<MultiVector> X   = MultiVectorFactory::Build(map, 1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

      X->putScalar(1.0);
      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

      Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

      // Use AMG directly as an iterative method
      {
        X->putScalar((Scalar)0.0);

        H->Iterate(*RHS, *X, its);

        X->norm2(norms);
        if (comm->getRank() == 0)
          out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
        results[run] = norms[0];
      }
    }

    TEST_EQUALITY(results[0] - results[1] < 1e-10, true);  // check results of EPETRA vs TPETRA
  }                                                        // comm->getSize == 1
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2)." << std::endl;
#endif

}  // SaPFactory_EpetraVsTpetra

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, ConstrainRowOptimalScalarPDE, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Don't test for complex - matrix reader won't work
  if (STS::isComplex) {
    success = true;
    return;
  }

  RCP<Matrix> P;
  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
  P                         = Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/SaP_constrainTest_P.mat", lib, comm);

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  sapFactory->optimalSatisfyPConstraintsForScalarPDEs(P);

  // check that row sums are all one by checking the norm of the vector
  // note: optimalSatisfyPConstraintsForScalarPDEs preserves row sum of original P (one in this case),
  //       but SatisfyPConstraints normalizes each row sum to one.
  RCP<MultiVector> X     = MultiVectorFactory::Build(P->getDomainMap(), 1);
  RCP<MultiVector> Bfact = MultiVectorFactory::Build(P->getRangeMap(), 1);
  X->putScalar(one);
  P->apply(*X, *Bfact, Teuchos::NO_TRANS, one, zero);
  Array<typename STS::magnitudeType> norms(1);
  Bfact->norm2(norms);
  out << "|| B_factory ones || = " << norms[0] << std::endl;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;
  TEST_FLOATING_EQUALITY(STS::magnitude(norms[0]), STS::magnitude(6.0), 100 * TMT::eps());

  // check that the min and max of each row are in [0,1]
  bool lowerViolation = false;
  bool upperViolation = false;
  for (size_t i = 0; i < (size_t)(P->getRowMap()->getLocalNumElements()); i++) {
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC> vals;
    P->getLocalRowView((LO)i, indices, vals);
    size_t nnz = indices.size();
    for (size_t j = 0; j < nnz; j++) {
      if (Teuchos::ScalarTraits<SC>::real(vals[j]) < Teuchos::ScalarTraits<SC>::real(zero)) lowerViolation = true;
      if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(one)) upperViolation = true;
      // if (STS::magnitude(vals[j]) < STS::magnitude(zero)) lowerViolation = true;
      // if (STS::magnitude(vals[j]) > STS::magnitude(one))  upperViolation = true;
    }
  }
  TEST_EQUALITY(lowerViolation, false);
  TEST_EQUALITY(upperViolation, false);

}  // SaPFactory_ConstrainRowOptimalScalarPDE

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, Test0, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, EpetraVsTpetra, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, ConstrainRowOptimalScalarPDE, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
