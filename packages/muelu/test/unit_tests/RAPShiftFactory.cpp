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

#include "MueLu_Utilities.hpp"
#include "MueLu_RAPShiftFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"

namespace MueLuTests {


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPShiftFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<RAPShiftFactory> rapFactory = rcp(new RAPShiftFactory);
    TEST_EQUALITY(rapFactory != Teuchos::null, true);

    out << *rapFactory << std::endl;
  } // Constructor test

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPShiftFactory, Correctness, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel, coarseLevel; TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    GO nx = 27*comm->getSize();
    RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
    fineLevel.Set("A",Op);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", Op->getRowMap(), galeriList);
    fineLevel.Set("Coordinates", coordinates);

    // Set "K" and "M" to be copies of A
    fineLevel.Set("K",Op);
    fineLevel.Set("M",Op);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory;
    sapFactory.SetFactory("P",rcpFromRef(tentpFactory));
    TransPFactory transPFactory;
    transPFactory.SetFactory("P", rcpFromRef(sapFactory));  //todo:rcpFromRef

    coarseLevel.Request("P",&sapFactory);
    coarseLevel.Request("R",&transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel,coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);

    RAPShiftFactory rap;
    Teuchos::ParameterList rapList;
    rapList.set("rap: shift", 1.0);
    rap.SetParameterList(rapList);

    rap.SetFactory("P", rcpFromRef(sapFactory));
    rap.SetFactory("R", rcpFromRef(transPFactory));

    coarseLevel.Request(rap);

    coarseLevel.Request("A",&rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A");
    RCP<Matrix> P = coarseLevel.Get< RCP<Matrix> >("P", &sapFactory);
    RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", &transPFactory);

    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(R->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = 2*R*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(Scalar)2.0,(Scalar)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    R->apply(*workVec2,*result1,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);

    RCP<Matrix> coarseOp = coarseLevel.Get< RCP<Matrix> >("A", &rap);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(R->getRangeMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);

    Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X))." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } // Correctness test

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPShiftFactory, ImplicitTranspose, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    typedef typename Teuchos::ScalarTraits<Scalar> TST;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    if (comm->getSize() > 1 && TestHelpers::Parameters::getLib() == Xpetra::UseEpetra ) {
      out << "Skipping ImplicitTranspose test for Epetra and #proc>1" << std::endl;
      return;
    }

    // build test-specific default factory manager
    RCP<FactoryManager> defManager = rcp(new FactoryManager());
    defManager->SetKokkosRefactor(false);
    defManager->SetFactory("A", rcp(MueLu::NoFactory::get(),false));         // dummy factory for A
    defManager->SetFactory("Nullspace", rcp(new NullspaceFactory()));        // real null space factory for Ptent
    defManager->SetFactory("Graph", rcp(new CoalesceDropFactory()));         // real graph factory for Ptent
    defManager->SetFactory("Aggregates", rcp(new CoupledAggregationFactory()));   // real aggregation factory for Ptent

    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    // overwrite default factory manager
    fineLevel.SetFactoryManager(defManager);
    coarseLevel.SetFactoryManager(defManager);

    GO nx = 19*comm->getSize();
    RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
    fineLevel.Set("A",Op);
    fineLevel.Set("K",Op);
    fineLevel.Set("M",Op);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", Op->getRowMap(), galeriList);
    fineLevel.Set("Coordinates", coordinates);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory;
    sapFactory.SetFactory("P",rcpFromRef(tentpFactory));
    TransPFactory transPFactory;
    transPFactory.SetFactory("P", rcpFromRef(sapFactory));
    coarseLevel.Request("P", &sapFactory);
    coarseLevel.Request("R", &transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel, coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);
    RAPShiftFactory rap;

    Teuchos::ParameterList rapList = *(rap.GetValidParameterList());
    rapList.set("transpose: use implicit", true);
    rapList.set("rap: shift",1.0);
    rap.SetParameterList(rapList);
    rap.SetFactory("P", rcpFromRef(sapFactory));
    rap.SetFactory("R", rcpFromRef(transPFactory));
    coarseLevel.Request("A", &rap);

    coarseLevel.Request(rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A");
    RCP<Matrix> P = coarseLevel.Get< RCP<Matrix> >("P", &sapFactory);
    RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", &transPFactory);


    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = 2*P^T*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(Scalar)2.0,(Scalar)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    P->apply(*workVec2,*result1,Teuchos::TRANS,(Scalar)1.0,(Scalar)0.0);

    RCP<Matrix> coarseOp = coarseLevel.Get< RCP<Matrix> >("A", &rap);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);

    Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the (Non-)Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit transpose of P." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } // ImplicitTranspose test



/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPShiftFactory,CreatePreconditioner_Factory, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
    typedef Node  NO;
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    using Teuchos::RCP;

    typedef typename Teuchos::ScalarTraits<Scalar> TST;
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    GO nx = 1999*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // The ugly way of managing input decks
    std::string myInputString =
"<ParameterList name=\"MueLu\">"
"  <ParameterList name=\"Factories\">"
"   <ParameterList name=\"myFilteredAFact\">"
"     <Parameter name=\"factory\"                             type=\"string\" value=\"FilteredAFactory\"/>"
"   </ParameterList>"
"   <ParameterList name=\"myProlongatorFact\">"
"     <Parameter name=\"factory\"                             type=\"string\" value=\"SaPFactory\"/>"
"     <Parameter name=\"A\"                                   type=\"string\" value=\"myFilteredAFact\"/>"
"   </ParameterList>"
"   <ParameterList name=\"myRAPFact\">"
"     <Parameter name=\"factory\"                             type=\"string\" value=\"RAPShiftFactory\"/>"
"     <Parameter name=\"P\"                                   type=\"string\" value=\"myProlongatorFact\"/>"
"     <Parameter name=\"rap: shift\"                          type=\"double\" value=\"1.0\"/>"
"   </ParameterList>"
"   <ParameterList name=\"mySmoo\">"
"     <Parameter name=\"factory\"                             type=\"string\" value=\"TrilinosSmoother\"/>"
"     <Parameter name=\"type\"                                type=\"string\" value=\"RELAXATION\"/>"
"   </ParameterList>"
" </ParameterList>"
" <ParameterList name=\"Hierarchy\">"
"   <Parameter name=\"max levels\"                            type=\"int\"      value=\"10\"/>"
"   <Parameter name=\"coarse: max size\"                      type=\"int\"      value=\"1000\"/>"
"   <Parameter name=\"verbosity\"                             type=\"string\"   value=\"None\"/>"
"   <Parameter name=\"use kokkos refactor\"                   type=\"bool\"     value=\"false\"/>"
"   <ParameterList name=\"All\">"
"     <Parameter name=\"startLevel\"                          type=\"int\"      value=\"0\"/>"
"     <Parameter name=\"Smoother\"                            type=\"string\"   value=\"mySmoo\"/>"
"     <Parameter name=\"A\"                                   type=\"string\"   value=\"myRAPFact\"/>"
"     <Parameter name=\"P\"                                   type=\"string\"   value=\"myProlongatorFact\"/>"
"     <Parameter name=\"CoarseSolver\"                        type=\"string\"   value=\"DirectSolver\"/>"
"   </ParameterList>"
" </ParameterList>"
"</ParameterList>"
     ;

    // Get the actual parameter list
    RCP<Teuchos::ParameterList> Params = Teuchos::getParametersFromXmlString(myInputString);
    Params->set("verbosity","none");

    Teuchos::ParameterList & pLevel0 = Params->sublist("level 0");
    pLevel0.set("K",A);
    pLevel0.set("M",A);

    // Build hierarchy
    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,*Params, coordinates, Teuchos::null);

    // Ready the test vector
    RCP<Level> Level1 = H->GetLevel(1);
    RCP<Matrix> P = Level1->Get< RCP<Matrix> >("P");
    RCP<Matrix> R = Level1->Get< RCP<Matrix> >("R");
    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(A->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = 2*P^T*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(Scalar)2.0,(Scalar)0.0);
    A->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    P->apply(*workVec2,*result1,Teuchos::TRANS,(Scalar)1.0,(Scalar)0.0);

    RCP<Matrix> coarseOp = Level1->Get< RCP<Matrix> >("A");

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);

    Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the (Non-)Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit transpose of P." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPShiftFactory,CreatePreconditioner_Easy, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
    typedef Node  NO;
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    using Teuchos::RCP;

    typedef typename Teuchos::ScalarTraits<Scalar> TST;
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    GO nx = 1999*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

    // The pretty way of managing input decks
    RCP<Teuchos::ParameterList> Params = rcp(new Teuchos::ParameterList);
    Params->set("rap: algorithm","shift");
    Params->set("rap: shift",1.0);
    Params->set("coarse: max size",1000);
    Params->set("verbosity","none");
    Params->set("use kokkos refactor",false);

    Teuchos::ParameterList & pLevel0 = Params->sublist("level 0");
    pLevel0.set("K",A);
    pLevel0.set("M",A);

    // Build hierarchy
    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,*Params,coordinates,Teuchos::null);

    // Ready the test vector
    RCP<Level> Level1 = H->GetLevel(1);
    RCP<Matrix> P = Level1->Get< RCP<Matrix> >("P");
    RCP<Matrix> R = Level1->Get< RCP<Matrix> >("R");
    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(A->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = 2*P^T*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(Scalar)2.0,(Scalar)0.0);
    A->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
    P->apply(*workVec2,*result1,Teuchos::TRANS,(Scalar)1.0,(Scalar)0.0);

    RCP<Matrix> coarseOp = Level1->Get< RCP<Matrix> >("A");

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);

    Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the (Non-)Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit transpose of P." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  }


#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPShiftFactory,Constructor,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPShiftFactory,Correctness,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPShiftFactory,ImplicitTranspose,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPShiftFactory,CreatePreconditioner_Factory,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPShiftFactory,CreatePreconditioner_Easy,Scalar,LO,GO,Node)
#include <MueLu_ETI_4arg.hpp>

} // namespace MueLuTests

