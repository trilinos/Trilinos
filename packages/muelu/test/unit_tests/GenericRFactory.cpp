/*
 * GenericRFactory.cpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */


#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

#include "MueLu_SaPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace MueLuTests {


  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(GenericRFactory, Constructor)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    out << "version: " << MueLu::Version() << std::endl;

    //TEST_EQUALITY(needs != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST(GenericRFactory, SymmetricProblem)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // generate problem
    LO maxLevels = 3;
    LO its=10;
    LO nEle = 63;
    const RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);
    RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>("Laplace1D", map, matrixParameters);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);
    Teuchos::Array<ST::magnitudeType> norms(1);
    nullSpace->norm1(norms);
    if (comm->getRank() == 0)
      out << "||NS|| = " << norms[0] << std::endl;

    // setup finest level
    RCP<Level> Finest = rcp( new Level() );
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

    // prepare default factory handler for graph
    RCP<DefaultFactoryHandlerBase> defHandler = rcp(new DefaultFactoryHandlerBase());
    defHandler->SetDefaultFactory("Graph",rcp(new CoalesceDropFactory()));

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy(defHandler) );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    H->SetLevel(Finest); // first associate level with hierarchy (for defaultFactoryHandler!)

    // define transfer operators
    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory(UCAggFact));
    RCP<SaPFactory>         Pfact = rcp( new SaPFactory(Ptentfact));
    RCP<RFactory>           Rfact = rcp( new GenericRFactory(Pfact) );
    RCP<GenericPRFactory>  PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    PRfact->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(TestHelpers::Parameters::getLib(), "RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::ParameterList status;
    status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);

    SmootherFactory coarseSolveFact(smooProto);
    H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

    H->GetLevel(1)->print(std::cout);
    H->GetLevel(2)->print(std::cout);
    H->GetLevel(3)->print(std::cout);

    // test some basic multgrid data
    RCP<Level> coarseLevel = H->GetLevel(2);
    RCP<Operator> P1 = coarseLevel->Get< RCP<Operator> >("P",NULL);
    RCP<Operator> R1 = coarseLevel->Get< RCP<Operator> >("R",NULL);
    RCP<Level> coarseLevel2 = H->GetLevel(3);
    RCP<Operator> P2 = coarseLevel2->Get< RCP<Operator> >("P",NULL);
    RCP<Operator> R2 = coarseLevel2->Get< RCP<Operator> >("R",NULL);

    Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > RP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(R1,false,P1,false);

    RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    RP->describe(*fos,Teuchos::VERB_EXTREME);

    //Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P1,true,P1,false);
    //PtP->describe(*fos,Teuchos::VERB_EXTREME);



    //R1->describe(*fos,Teuchos::VERB_EXTREME);

    /*RCP<CrsOperator> crsP1 = rcp_dynamic_cast<CrsOperator>(P1);
    RCP<CrsMatrix> crsMat = crsP1->getCrsMatrix();
    RCP<Xpetra::EpetraCrsMatrix> epcrsMat = rcp_dynamic_cast<Xpetra::EpetraCrsMatrix>(crsMat);
    RCP<Epetra_CrsMatrix> epMat = epcrsMat->getEpetra_CrsMatrixNonConst();
    EpetraExt::RowMatrixToMatrixMarketFile( "Test.mat", *epMat);*/

    //P1->describe(*fos,Teuchos::VERB_EXTREME);

    //R1->getRangeMap()->describe(*fos,Teuchos::VERB_EXTREME);
    //P1->describe(*fos,Teuchos::VERB_EXTREME);
    //R1->describe(*fos,Teuchos::VERB_EXTREME);

    // todo test me
    /*RCP<Operator> R1T = MueLu::Utils2<SC,LO,GO>::Transpose(P1,true);
    RCP<Operator> R2T = MueLu::Utils2<SC,LO,GO>::Transpose(P2,true);*/

    /*RCP<Vector> X1 = VectorFactory::Build(P1->getDomainMap());
    RCP<Vector> X2 = VectorFactory::Build(R1->getRangeMap());
    RCP<Vector> Y1 = VectorFactory::Build(P1->getRangeMap());
    RCP<Vector> Y2 = VectorFactory::Build(R1->getDomainMap());

    X1->putScalar(1.0);
    X2->putScalar(1.0);

    P1->apply(*X1,*Y1);
    R1T->apply(*X2,*Y2);
    std::cout << Y1->norm1() << std::endl;
    std::cout << Y2->norm1() << std::endl;

    std::cout << X1->norm1() << std::endl;
    std::cout << X2->norm1() << std::endl;*/

    //Y1->update(-1.0, *Y2, 1.0);
    //std::cout << Y1->norm1() << std::endl;
    //std::cout << Y1->norm2() << std::endl;
    //TEST_EQUALITY(Y1->norm1() < 1e-6, true);

  }


}//namespace MueLuTests


