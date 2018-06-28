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
// // Questions? Contact
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

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_AmesosSmoother.hpp>
#include <MueLu_AmesosSmoother.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_FactoryManagerBase.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_PFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_KOKKOSCORE
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#endif

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Hierarchy> H = rcp(new Hierarchy);

    TEST_INEQUALITY(H, Teuchos::null);

  } //Constructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, DescriptionCaching, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    /*
     * Test to confirm that Hierarchy::description():
     * a) gives the same result before and after a call to SetupRe() when the number of levels has not changed
     * b) gives a different result before and after a call to SetupRe() when the number of levels has changed
     *
     * The occasion for this test is the introduction of caching for the result of description(), to avoid redundant
     * computation.
     */
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int numRows = 399;
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(numRows);
    GO nx = numRows;
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    Teuchos::ParameterList MueLuList;
    MueLuList.set("verbosity","none");
    //MueLuList.set("verbosity","high");
    MueLuList.set("coarse: max size",numRows-1); // make it so we want two levels
    MueLuList.set("max levels",2);

    Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > H =
        MueLu::CreateXpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
            A, MueLuList, coordinates, Teuchos::null);

    // confirm that we did get a hierarchy with two levels -- a sanity check for this test
    TEST_EQUALITY(2, H->GetGlobalNumLevels());

    using namespace std;
    string descriptionTwoLevel = H->description();

    // SetupRe() will reset description, but since we haven't changed anything, the result from description() should not change
    H->SetupRe();
    string descriptionActual = H->description();
    TEST_EQUALITY(descriptionActual,descriptionTwoLevel);

    // now, we allow a larger coarse size; we should get just one level during SetupRe()
    H->SetMaxCoarseSize(numRows+1);
    H->SetupRe();
    // as a sanity check for the test, make sure that we do have just one level
    TEST_EQUALITY(1, H->GetGlobalNumLevels());
    descriptionActual = H->description();

    // since the number of levels has changed, the new description should differ from the old
    TEST_INEQUALITY(descriptionActual, descriptionTwoLevel);
  }//DescriptionCaching

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetAndGetLevel, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    Hierarchy H;
    //   RCP<Level> level = rcp(new Level());
    //   H.SetLevel(level);
    //   RCP<Level> dupLevel = H.GetLevel(1);

    //  TEST_EQUALITY(level, dupLevel);

  }//SetAndGetLevel

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, GetNumLevels, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    Hierarchy H;
    RCP<Level> levelOne = rcp(new Level() );
    RCP<Level> levelTwo = rcp(new Level() );
    RCP<Level> levelThree = rcp(new Level() );
    //   H.SetLevel(levelOne);
    //   H.SetLevel(levelTwo);
    //   H.SetLevel(levelThree);
    //   TEST_EQUALITY(H.GetNumLevels(), 3);

  }//GetNumLevels

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, KeepAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 399*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    Hierarchy H(A);
    H.SetMaxCoarseSize(1);
    H.GetLevel(0)->Set("Coordinates", coordinates);

    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", Teuchos::null);
    M.SetFactory("CoarseSolver", Teuchos::null);

    H.GetLevel(0)->Keep("Aggregates", CoupledAggFact.get());
    H.Setup(M, 0, 2);

    for (LocalOrdinal l=0; l<H.GetNumLevels()-1;l++) {
      TEST_EQUALITY(H.GetLevel(l)->IsAvailable("Aggregates", CoupledAggFact.get()), true);
    }

  } //FullPopulate_KeepAggregates

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, Iterate, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<Scalar> TST;

#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    out << "version: " << MueLu::Version() << std::endl;

    //matrix
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 6561*comm->getSize();  //=8*3^6
    RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);
    RCP<const Map > map = Op->getRowMap();

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", map, galeriList);

    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
    nullSpace->putScalar( (Scalar) 1.0);
    Teuchos::Array<typename TST::magnitudeType> norms(1);
    nullSpace->norm1(norms);

    MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> H;
    H.setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<MueLu::Level> Finest = H.GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    Finest->Set("Coordinates", coordinates);
    Finest->Set("Nullspace", nullSpace);
    Finest->Set("A", Op);

    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<CoalesceDropFactory> cdFact;
    RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory());

    RCP<SaPFactory>    Pfact  = rcp( new SaPFactory() );
    RCP<TransPFactory> Rfact  = rcp( new TransPFactory());
    RCP<RAPFactory>    Acfact = rcp( new RAPFactory() );

    RCP<SmootherPrototype> smooProto = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSmootherPrototype("Gauss-Seidel", 2);
    RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherPrototype> coarseProto = rcp(new DirectSolver("Klu"));
    RCP<SmootherFactory> coarseSolveFact = rcp( new SmootherFactory(coarseProto, Teuchos::null));

    int maxLevels = 5;

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", TentPFact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);
    M.SetFactory("Coordinates", TentPFact);

    H.Setup(M, 0, maxLevels);

    RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

    X->setSeed(846930886);
    X->randomize();
    //Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

    X->norm2(norms);
    X->scale(1/norms[0]);
    X->norm2(norms);
    out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    RHS->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);

    X->norm2(norms);
    out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(std::ios::fixed) <<
      std::setprecision(10) << norms[0] << std::endl;

    norms = Utilities::ResidualNorm(*Op, *X, *RHS);
    out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
    TEST_EQUALITY(norms[0]<1e-10, true);

  } //Iterate

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, IterateWithImplicitRestriction, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    //matrix
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 6561*comm->getSize();  //=8*3^6
    RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);
    RCP<const Map > map = Op->getRowMap();

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", map, galeriList);

    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
    nullSpace->putScalar( (Scalar) 1.0);
    Teuchos::Array<typename TST::magnitudeType> norms(1);
    nullSpace->norm1(norms);

    MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> H;
    H.SetImplicitTranspose(true);
    H.setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<MueLu::Level> Finest = H.GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Op);
    Finest->Set("Nullspace", nullSpace);
    Finest->Set("Coordinates", coordinates);

    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);
    RCP<CoalesceDropFactory> cdFact;
    RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory());

    RCP<SaPFactory>         Pfact = rcp( new SaPFactory() );
    RCP<TransPFactory>      Rfact = rcp( new TransPFactory());
    RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
    Teuchos::ParameterList Aclist = *(Acfact->GetValidParameterList());
    Aclist.set("transpose: use implicit", true);
    Acfact->SetParameterList(Aclist);

    RCP<SmootherPrototype> smooProto = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSmootherPrototype("Gauss-Seidel", 2);
    RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherPrototype> coarseProto = rcp(new DirectSolver("Klu"));
    RCP<SmootherFactory> coarseSolveFact = rcp( new SmootherFactory(coarseProto, Teuchos::null));

    int maxLevels = 5;

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", TentPFact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);
    M.SetFactory("Coordinates", TentPFact);

    H.Setup(M, 0, maxLevels);

    //  Teuchos::ParameterList status;
    // status.print(out, Teuchos::ParameterList::PrintOptions().indent(2));

    RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

    X->setSeed(846930886);
    X->randomize();
    //Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

    X->norm2(norms);
    X->scale(1/norms[0]);
    X->norm2(norms);
    out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    RHS->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);

    X->norm2(norms);
    out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(std::ios::fixed) <<
      std::setprecision(10) << norms[0] << std::endl;

    norms = Utilities::ResidualNorm(*Op, *X, *RHS);
    out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
    TEST_EQUALITY(norms[0]<1e-10, true);

  } //Iterate

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchy1level, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(299*comm->getSize());

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.setVerbLevel(Teuchos::VERB_HIGH);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);

    bool r = H.Setup(0, Teuchos::null, rcpFromRef(M0), Teuchos::null); TEST_EQUALITY(r, true); // cf. Teuchos Bug 5214

    RCP<Level> l0 = H.GetLevel(0);

    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    //TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchy1levelv2, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    //MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra);
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    MUELU_TESTING_SET_OSTREAM;
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(299*comm->getSize());

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.setVerbLevel(Teuchos::VERB_HIGH);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);

    H.Setup(M0, 0, 1);

    RCP<Level> l0 = H.GetLevel(0);

    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    //TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations = 10;
    H.Iterate(*RHS, *X, iterations);
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchy2level, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra);
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(299*comm->getSize());

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.setVerbLevel(Teuchos::VERB_HIGH);
    H.SetMaxCoarseSize(50);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);

    FactoryManager M1; // first coarse level (Plain aggregation)
    M1.SetKokkosRefactor(false);
    M1.SetFactory("A", rcp(new RAPFactory()));
    M1.SetFactory("P", rcp(new TentativePFactory()));

    FactoryManager M2; // last level (SA)
    M2.SetKokkosRefactor(false);
    M2.SetFactory("A", rcp(new RAPFactory()));
    M2.SetFactory("P", rcp(new SaPFactory()));

    bool r; // cf. Teuchos Bug 5214
    r = H.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1));  TEST_EQUALITY(r, false);
    r = H.Setup(1, rcpFromRef(M0), rcpFromRef(M1), Teuchos::null); TEST_EQUALITY(r, true);

    RCP<Level> l0 = H.GetLevel(0);
    RCP<Level> l1 = H.GetLevel(1);

    /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      l0->print(*stdout,Teuchos::VERB_EXTREME);
      l1->print(*stdout,Teuchos::VERB_EXTREME);*/

    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

    TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    // TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);
#  endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchy3level, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    MUELU_TESTING_SET_OSTREAM;

#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 299*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.setVerbLevel(Teuchos::VERB_HIGH);
    H.SetMaxCoarseSize(50);
    H.GetLevel(0)->Set("Coordinates", coordinates);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);

    FactoryManager M1; // first coarse level (Plain aggregation)
    M1.SetKokkosRefactor(false);
    M1.SetFactory("A", rcp(new RAPFactory()));
    RCP<FactoryBase> P = rcp(new TentativePFactory());
    M1.SetFactory("P", P);
    M1.SetFactory("Ptent", P); //FIXME: can it be done automatically in FactoryManager?

    FactoryManager M2; // last level (SA)
    M2.SetKokkosRefactor(false);
    M2.SetFactory("A", rcp(new RAPFactory()));
    M2.SetFactory("P", rcp(new SaPFactory()));

    bool r; // cf. bug Teuchos Bug 5214
    r = H.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1)); TEST_EQUALITY(r, false);
    r = H.Setup(1, rcpFromRef(M0), rcpFromRef(M1), rcpFromRef(M2));   TEST_EQUALITY(r, false);
    r = H.Setup(2, rcpFromRef(M1), rcpFromRef(M2), Teuchos::null ); TEST_EQUALITY(r, true);

    RCP<Level> l0 = H.GetLevel(0);
    RCP<Level> l1 = H.GetLevel(1);
    RCP<Level> l2 = H.GetLevel(2);

    /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      l0->print(*stdout,Teuchos::VERB_EXTREME);
      l1->print(*stdout,Teuchos::VERB_EXTREME);
      l2->print(*stdout,Teuchos::VERB_EXTREME);*/

    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

    TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

    TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    // TEST_EQUALITY(l2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchy3levelFacManagers, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    //MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 299*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.setVerbLevel(Teuchos::VERB_HIGH);
    H.SetMaxCoarseSize(50);
    H.GetLevel(0)->Set("Coordinates", coordinates);

    // setup smoother factory
    RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LocalOrdinal) 2);
    ifpackList.set("relaxation: damping factor", (Scalar) 0.9); // 0.7
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = Teuchos::rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, MueLu::Exceptions::Incompatible, "MueLu: UnitTest::Hierarchy::SetupHierarchy3levelFacManagers: Dynamic cast to IfpackSmoother failed. Epetra needs Scalar=double, LocalOrdinal=GlobalOrdinal=int and Node=Serial.");

    RCP<SmootherFactory> preSmooFact;
    RCP<SmootherFactory> postSmooFact;
    preSmooFact = rcp( new SmootherFactory(smooProto) );
    postSmooFact = rcp( new SmootherFactory(smooProto) );
    preSmooFact->SetSmootherPrototypes(smooProto,Teuchos::null );
    postSmooFact->SetSmootherPrototypes(Teuchos::null,smooProto);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);
    M0.SetFactory("Smoother", preSmooFact);

    FactoryManager M1; // first coarse level (Plain aggregation)
    M1.SetKokkosRefactor(false);
    M1.SetFactory("A", rcp(new RAPFactory()));
    RCP<FactoryBase> PFact = rcp(new TentativePFactory());
    M1.SetFactory("P", PFact);
    M1.SetFactory("Ptent", PFact); //FIXME: can it be done automatically in FactoryManager?
    M1.SetFactory("Smoother", postSmooFact);

    FactoryManager M2; // last level (SA)
    M2.SetKokkosRefactor(false);
    M2.SetFactory("A", rcp(new RAPFactory()));
    M2.SetFactory("P", rcp(new SaPFactory()));

    bool r; // cf. bug Teuchos Bug 5214
    H.EnableGraphDumping("hierarchy_test_graph",0);

    r = H.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1)); TEST_EQUALITY(r, false);
    r = H.Setup(1, rcpFromRef(M0), rcpFromRef(M1), rcpFromRef(M2));   TEST_EQUALITY(r, false);
    r = H.Setup(2, rcpFromRef(M1), rcpFromRef(M2), Teuchos::null);  TEST_EQUALITY(r, true);

    RCP<Level> l0 = H.GetLevel(0);
    RCP<Level> l1 = H.GetLevel(1);
    RCP<Level> l2 = H.GetLevel(2);

    /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      l0->print(*stdout,Teuchos::VERB_EXTREME);
      l1->print(*stdout,Teuchos::VERB_EXTREME);
      l2->print(*stdout,Teuchos::VERB_EXTREME);*/

    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);

    TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

    TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    // TEST_EQUALITY(l2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, SetupHierarchyTestBreakCondition, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx  = 299*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.SetMaxCoarseSize(299*comm->getSize()); // set max coarse size to fit problem size (-> 1 level method)
    H.setVerbLevel(Teuchos::VERB_HIGH);
    H.GetLevel(0)->Set("Coordinates", coordinates);

    // Multigrid setup phase (using default parameters)
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);

    FactoryManager M1; // first coarse level (Plain aggregation)
    M1.SetKokkosRefactor(false);
    M1.SetFactory("A", rcp(new RAPFactory()));
    M1.SetFactory("P", rcp(new TentativePFactory()));

    bool r; // cf. bug Teuchos Bug 5214
    r = H.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1)); TEST_EQUALITY(r, true);
    TEST_EQUALITY(H.GetNumLevels(),1);

    RCP<Level> l0 = H.GetLevel(0);
    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
    TEST_EQUALITY(l0->IsAvailable("P",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("R",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    // TEST_EQUALITY(l0->GetKeepFlag("PostSmoother",  MueLu::NoFactory::get()), MueLu::Final); //direct solve

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();

    X->putScalar( (Scalar) 0.0);

    int iterations=10;
    H.Iterate(*RHS, *X, iterations);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, Write, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx = 30;
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(nx);
    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // Multigrid Hierarchy
    Hierarchy H(A);
    H.SetDefaultVerbLevel(MueLu::Low);
    H.SetMaxCoarseSize(29);
    H.GetLevel(0)->Set("Coordinates", coordinates);

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("Smoother", Teuchos::null);
    M.SetFactory("CoarseSolver", Teuchos::null);
    H.Setup(M, 0, 2);

    TEST_THROW( H.Write(1,0), MueLu::Exceptions::RuntimeError );    //start level is greater than end level
    TEST_THROW( H.Write(0,1000), MueLu::Exceptions::RuntimeError ); //end level is too big

    // Write matrices out, read fine A back in, and check that the read was ok
    // by using a matvec with a random vector.
    // JJH: 22-Feb-2016 Append scalar type to file name. The theory is that for dashboard
    //      tests with multiple Scalar instantiations of this test, a test with Scalar type
    //      A could try to read in the results of the test with Scalar type B, simply because
    //      the test with type B overwrote A's output matrix file.  A better solution would be
    //      to write to a file stream, but this would involve writing new interfaces to Epetra's
    //      file I/O capabilities.
    std::string tname = typeid(Scalar).name();
    tname = tname + typeid(LocalOrdinal).name();
    tname = tname + typeid(GlobalOrdinal).name();
#ifdef HAVE_MUELU_KOKKOSCORE
    std::string nn = Kokkos::Compat::KokkosDeviceWrapperNode<typename Node::execution_space>::name();
    nn.erase(std::remove(nn.begin(), nn.end(), '/'), nn.end());
    tname = tname + nn;
#endif
    tname = "_" + tname;
    LocalOrdinal zero = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
    //Only write out the fine level matrix, since that is the only data file we test against.
    H.Write(zero,zero,tname);

    Teuchos::Array<typename TST::magnitudeType> norms(1);

    out << "random status: " << rand() << std::endl;
    std::srand(595343843);
    std::rand();
    std::string infile = "A_0" + tname + ".m";
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<Matrix> Ain = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(infile, lib, comm);
    remove(infile.c_str());
    infile = "colmap_A_0" + tname + ".m";    remove(infile.c_str());
    infile = "domainmap_A_0" + tname + ".m"; remove(infile.c_str());
    infile = "rangemap_A_0" + tname + ".m";  remove(infile.c_str());
    infile = "rowmap_A_0" + tname + ".m";    remove(infile.c_str());
    RCP<Vector> randomVec = VectorFactory::Build(A->getDomainMap(),false);
    randomVec->randomize();
    out << "randomVec norm: " << randomVec->norm2() << std::endl;
    RCP<Vector> A_v = VectorFactory::Build(A->getRangeMap(),false);
    A->apply(*randomVec,*A_v,Teuchos::NO_TRANS,1,0);
    out << "A_v norm: " << A_v->norm2() << std::endl;

    RCP<Vector> Ain_v = VectorFactory::Build(Ain->getRangeMap(),false);
    Ain->apply(*randomVec,*Ain_v,Teuchos::NO_TRANS,1,0);
    out << "Ain_v norm: " << Ain_v->norm2() << std::endl;

    RCP<MultiVector> diff = VectorFactory::Build(A->getRangeMap());
    //diff = A_v + (-1.0)*(Ain_v) + 0*diff
    diff->update(1.0,*A_v,-1.0,*Ain_v,0.0);

    diff->norm2(norms);
    out << "||diff|| = " << norms[0] << std::endl;
    TEST_EQUALITY(norms[0]<1e-15, true);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Hierarchy, BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra);

    out << "===== Generating matrices =====" << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Teuchos::ParameterList matrixList;
    matrixList.set<GO>("nx",          100);
    matrixList.set<GO>("ny",          100);
    matrixList.set("matrixType",  "Laplace2D");

    // Construct block matrix
    RCP<Matrix> A = TestHelpers::TpetraTestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBlockMatrix(matrixList,Xpetra::UseTpetra);
    if(A==Teuchos::null) { // if A is Teuchos::null, we could not build the matrix as it is not instantiated in Tpetra
      out << "Skipping test" << std::endl;
      return;
    }

    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("2D", A->getRowMap(), matrixList);

    // extract information
    RCP<const Map>        rangeMap = A->getRangeMap();
    Xpetra::UnderlyingLib lib      = rangeMap->lib();

    // Construct our own point operators for level 1
    matrixList.set("nx", matrixList.get<GO>("nx")*3);
    matrixList.set("ny", matrixList.get<GO>("ny"));
    RCP<Matrix> A_point = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMatrix(matrixList, lib);;

    out << "===== Setting up helper hierarchy =====" << std::endl;

    // Setup helper hierarchy
    Hierarchy Haux(A_point);
    Haux.SetMaxCoarseSize(1);
    Haux.SetDefaultVerbLevel(MueLu::None);
    Haux.GetLevel(0)->Set("Coordinates", coordinates);

    FactoryManager Maux;
    Maux.SetKokkosRefactor(false);
    const FactoryBase* nullFactory  = Maux.GetFactory("Nullspace").get();

    Haux.Keep("Nullspace", nullFactory);

    // Level 1 nullspace is only constructed/requested by level 2 factories. So
    // we need at least three level hierarchy here, unless we somehow call
    // NullspaceFactory::Build after hierarchy setup
    Haux.Setup(Maux, 0, 3);
    TEST_EQUALITY(Haux.GetNumLevels(), 3);

    // Build block SGS smoother
    std::string ifpack2Type;
    Teuchos::ParameterList ifpack2List;
    ifpack2Type = "RBILUK";
    out << ifpack2Type << std::endl;
    RCP<SmootherPrototype> smooProto = Teuchos::rcp( new Ifpack2Smoother(ifpack2Type) );
    RCP<SmootherFactory>   SmooFact  = rcp( new SmootherFactory(smooProto) );

    // Setup hierarchy managers
    // We need three managers due to restriction in HierarchyHelpers: they
    // request R, A, and P for coarse level even if User provides them The quick
    // fix proposed here is to simply set corresponding fatories to NoFactory, so
    // they would fetch user data.
    FactoryManager M0, M1, M2;
    M0.SetKokkosRefactor(false);
    M0.SetFactory("Smoother", SmooFact);
    M1.SetKokkosRefactor(false);
    M1.SetFactory("A",        MueLu::NoFactory::getRCP());
    M1.SetFactory("P",        MueLu::NoFactory::getRCP());
    M1.SetFactory("R",        MueLu::NoFactory::getRCP());
    M2.SetKokkosRefactor(false);

    out << "===== Setting up mixed hierarchy =====" << std::endl;

    // Setup mixed hierarchy
    HierarchyManager mueluManager;
    mueluManager.AddFactoryManager(0, 1, Teuchos::rcpFromRef(M0));
    mueluManager.AddFactoryManager(1, 1, Teuchos::rcpFromRef(M1));
    mueluManager.AddFactoryManager(2, 1, Teuchos::rcpFromRef(M2));

    RCP<Hierarchy> Hrcp = mueluManager.CreateHierarchy();
    Hierarchy&     H    = *Hrcp;
    H.SetDefaultVerbLevel(MueLu::Low | MueLu::Debug);

    RCP<Level> l0 = H.GetLevel(0);
    l0->Set("A",           A);
    H.AddNewLevel();
    RCP<Level> l1     = H   .GetLevel(1);
    RCP<Level> l1_aux = Haux.GetLevel(1);
    l1->Set("P",          l1_aux->Get<RCP<Operator> >   ("P"));
    l1->Set("R",          l1_aux->Get<RCP<Operator> >   ("R"));
    l1->Set("A",          l1_aux->Get<RCP<Matrix> >     ("A"));
    l1->Set("Nullspace",  l1_aux->Get<RCP<MultiVector> >("Nullspace", nullFactory));

    mueluManager.SetupHierarchy(H);

    out << "===== Solving =====" << std::endl;

    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();
    X->putScalar( Teuchos::ScalarTraits<Scalar>::zero());

    int iterations = 10;
    H.IsPreconditioner(false);
    H.Iterate(*RHS, *X, iterations);
#endif
    TEST_EQUALITY(0,0);
  }


# define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, Constructor, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, DescriptionCaching, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetAndGetLevel, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, GetNumLevels, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, KeepAggregates, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, Iterate, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, IterateWithImplicitRestriction, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchy1level, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchy1levelv2, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchy2level, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchy3level, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchy3levelFacManagers, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, SetupHierarchyTestBreakCondition, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, Write, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Hierarchy, BlockCrs, Scalar, LO, GO, Node)

# include <MueLu_ETI_4arg.hpp>

}//namespace MueLuTests

//Note from JG:
// For UnitTest,  TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Matrix> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel");
// should be replaced by
// TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Matrix> >("PreSmoother"), preSmoother);
// testing if preSmoother->GetType() == "Ifpack: Gauss-Seidel" should be a unit test of the class IfpackSmoother


//TODO unit test:
// test if Hierarchy::Iterate(X,Y) works when X == Y (ie: do we need to test explicitly if X==Y and make a temporary copy inside of Iterate() ?)
