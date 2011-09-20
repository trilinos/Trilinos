#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_SaPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(SaPFactory, Test0)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    TEST_EQUALITY(sapFactory != Teuchos::null, true);

    out << *sapFactory << std::endl;

  }

  TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    sapFactory->SetDampingFactor( (Scalar)4/3 );
    TEST_EQUALITY(((Scalar)4/3) == sapFactory->GetDampingFactor(), true);
    sapFactory->SetDiagonalView("roomWithAView");
    TEST_EQUALITY( sapFactory->GetDiagonalView(), "roomWithAView");

  } //GetSetMethods

  TEUCHOS_UNIT_TEST(SaPFactory, SaPFactory_EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Compare results of Epetra and Tpetra" << std::endl;
    out << "for 3 level AMG solver using smoothed aggregation with" << std::endl;
    out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::Array<ST::magnitudeType> results(2);

    // run test only on 1 proc
    if(comm->getSize() == 1)
    {
        Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

        // run Epetra and Tpetra test
        for (int run = 0; run < 2; run++)
        {
            if (run == 0) lib = Xpetra::UseEpetra;
            else lib = Xpetra::UseTpetra;

            // generate problem
            LO maxLevels = 3;
            LO its=10;
            LO nEle = 63;
            const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
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
            RCP<RFactory>           Rfact = rcp( new TransPFactory() );
            RCP<GenericPRFactory>  PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
            RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
            PRfact->SetMaxCoarseSize(1);

            // setup smoothers
            Teuchos::ParameterList smootherParamList;
            smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
            smootherParamList.set("relaxation: sweeps", (LO) 1);
            smootherParamList.set("relaxation: damping factor", (SC) 1.0);
            RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(lib, "RELAXATION", smootherParamList) );
            RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
            Acfact->setVerbLevel(Teuchos::VERB_HIGH);

            Teuchos::ParameterList status;
            status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);

            SmootherFactory coarseSolveFact(smooProto);
            H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

            // test some basic multgrid data
            RCP<Level> coarseLevel = H->GetLevel(2);
            RCP<Operator> P1 = coarseLevel->Get< RCP<Operator> >("P",NULL);
            RCP<Operator> R1 = coarseLevel->Get< RCP<Operator> >("R",NULL);
            TEST_EQUALITY(P1->getGlobalNumRows(), 63);
            TEST_EQUALITY(P1->getGlobalNumCols(), 21);
            TEST_EQUALITY(R1->getGlobalNumRows(), 21);
            TEST_EQUALITY(R1->getGlobalNumCols(), 63);
            RCP<Level> coarseLevel2 = H->GetLevel(3);
            RCP<Operator> P2 = coarseLevel2->Get< RCP<Operator> >("P",NULL);
            RCP<Operator> R2 = coarseLevel2->Get< RCP<Operator> >("R",NULL);
            TEST_EQUALITY(P2->getGlobalNumRows(), 21);
            TEST_EQUALITY(P2->getGlobalNumCols(), 7);
            TEST_EQUALITY(R2->getGlobalNumRows(), 7);
            TEST_EQUALITY(R2->getGlobalNumCols(), 21);

            Teuchos::RCP<Xpetra::Operator<Scalar,LO,GO> > PtentTPtent = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P1,true,P1,false);
            TEST_EQUALITY(PtentTPtent->getGlobalMaxNumRowEntries()-3<1e-12, true);
            TEST_EQUALITY(P1->getGlobalMaxNumRowEntries()-2<1e-12, true);
            TEST_EQUALITY(P2->getGlobalMaxNumRowEntries()-2<1e-12, true);

            // Define RHS
            RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
            RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

            X->putScalar(1.0);
            X->norm2(norms);
            if (comm->getRank() == 0)
              out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

            Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

            // Use AMG directly as an iterative method
            {
              X->putScalar( (SC) 0.0);

              H->Iterate(*RHS,its,*X);

              X->norm2(norms);
              if (comm->getRank() == 0)
                out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
              results[run] = norms[0];
            }
        }

        TEST_EQUALITY(results[0] - results[1] < 1e-10, true); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1

  } //MakeTentative

}//namespace MueLuTests

