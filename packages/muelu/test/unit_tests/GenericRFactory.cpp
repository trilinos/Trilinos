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

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

    // define transfer operators
    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory(UCAggFact));
    RCP<SaPFactory>         Pfact = rcp( new SaPFactory(Ptentfact));
    RCP<RFactory>           Rfact = rcp( new GenericRFactory(Pfact) );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(TestHelpers::Parameters::getLib(), "RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::ParameterList status;
    status = H->FullPopulate(*Pfact,*Rfact, *Acfact,*SmooFact,0,maxLevels);

    SmootherFactory coarseSolveFact(smooProto);
    H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

    RCP<Level> coarseLevel = H->GetLevel(2);
    RCP<Operator> P1 = coarseLevel->Get< RCP<Operator> >("P");
    RCP<Operator> R1 = coarseLevel->Get< RCP<Operator> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(3);
    RCP<Operator> P2 = coarseLevel2->Get< RCP<Operator> >("P");
    RCP<Operator> R2 = coarseLevel2->Get< RCP<Operator> >("R");

    TEST_EQUALITY(Finest->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(Finest->IsAvailable("PostSmoother"), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother"), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother"), false);

    // test some basic multgrid data
    TEST_EQUALITY(P1->getGlobalNumEntries(), R1->getGlobalNumEntries());
    TEST_EQUALITY(P1->getGlobalNumRows(), R1->getGlobalNumCols());
    TEST_EQUALITY(P1->getGlobalNumCols(), R1->getGlobalNumRows());
    TEST_EQUALITY(P2->getGlobalNumEntries(), R2->getGlobalNumEntries());
    TEST_EQUALITY(P2->getGlobalNumRows(), R2->getGlobalNumCols());
    TEST_EQUALITY(P2->getGlobalNumCols(), R2->getGlobalNumRows());


    //RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    // since A is chosen symmetric, it is P^T = R
    // check P^T * P = R * P
    // note: the Epetra matrix-matrix multiplication using implicit transpose is buggy in parallel case
    //       (for multiplication of a square matrix with a rectangular matrix)
    //       however it seems to work for two rectangular matrices
    Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > RP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(R1,false,P1,false);
    Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P1,true,P1,false);

    RCP<Vector> x = VectorFactory::Build(RP->getDomainMap());
    RCP<Vector> bRP  = VectorFactory::Build(RP->getRangeMap());
    RCP<Vector> bPtP = VectorFactory::Build(PtP->getRangeMap());

    x->randomize();
    RP->apply(*x,*bRP);
    PtP->apply(*x,*bPtP);

    TEST_EQUALITY(bRP->norm1() - bPtP->norm1() < 1e-12, true);

    Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > RP2 = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(R2,false,P2,false);
    Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtP2 = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P2,true,P2,false);

    x = VectorFactory::Build(RP2->getDomainMap());
    bRP  = VectorFactory::Build(RP2->getRangeMap());
    bPtP = VectorFactory::Build(PtP2->getRangeMap());

    x->randomize();
    RP2->apply(*x,*bRP);
    PtP2->apply(*x,*bPtP);

    TEST_EQUALITY(bRP->norm1() - bPtP->norm1() < 1e-12, true);


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



  }


}//namespace MueLuTests


