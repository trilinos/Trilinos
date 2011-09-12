#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_PRFactory.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST(Hierarchy,Constructor)
{
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEST_INEQUALITY(H, Teuchos::null);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy,SetAndGetLevel)
{
  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> level = rcp(new Level());
  H.SetLevel(level);
  RCP<Level> dupLevel = H.GetLevel(1);

  TEST_EQUALITY(level, dupLevel);

}//SetAndGetLevel

TEUCHOS_UNIT_TEST(Hierarchy,NumberOfLevels)
{

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> levelOne = rcp(new Level() );
  RCP<Level> levelTwo = rcp(new Level() );
  RCP<Level> levelThree = rcp(new Level() );
  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);
  H.SetLevel(levelThree);
  TEST_EQUALITY(H.GetNumberOfLevels(), 3);

}//NumberOfLevels

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_NoFactoriesGiven)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  levelOne->Request("A");
  levelOne->Set("A",A);

  out << "Providing no factories to FillHierarchy." << std::endl;
  H.FillHierarchy();

  bool ToF = H.PrintResidualHistory();
  H.PrintResidualHistory(!ToF);
  TEST_INEQUALITY(H.PrintResidualHistory(), ToF);
} //FillHierarchy_NoFactoriesGiven

#ifdef FAILING_ONE_PROC_TERMINATING_BADLY_DONT_KNOW_WHY
TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PRFactoryOnly)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  H.SetLevel(levelOne);

  levelOne->Request("A");
  levelOne->Set("A",A);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just PR factory to FillHierarchy." << std::endl;
  Teuchos::ParameterList status;
  status = H.FillHierarchy(PRFact);
  TEST_EQUALITY(status.get("fine nnz",(Xpetra::global_size_t)-1), 295);
  TEST_EQUALITY(status.get("total nnz",(Xpetra::global_size_t)-1), 422);
  TEST_EQUALITY(status.get("start level",-1), 0);
  TEST_EQUALITY(status.get("end level",-1), 1);
  TEST_FLOATING_EQUALITY(status.get("operator complexity",(SC)-1.0),1.43051,1e-5);

} //FillHierarchy_PRFactoryOnly
#endif

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_BothFactories)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  H.SetLevel(levelOne);

  levelOne->Request("A");
  levelOne->Set("A",A);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<TransPFactory>  RFact = rcp(new TransPFactory(PFact));
  GenericPRFactory PRFact(PFact,RFact);
  RAPFactory    AcFact(rcpFromRef(PRFact));

  out << "Providing both factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact,AcFact);
} //FillHierarchy_BothFactories

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {


  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  RCP<Level> levelTwo = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);

  levelOne->Request("A");
  levelOne->Set("A",A);

//   H.SetSmoothers();

//   TEST_EQUALITY(H.GetLevel(1)->template Get< RCP<SmootherBase> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel");
//   TEST_EQUALITY(H.GetLevel(1)->template Get< RCP<SmootherBase> >("PostSmoother")->GetType(),"Ifpack: Gauss-Seidel");

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Jacobi");
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto) );
  H.SetSmoothers(*smooFactory);
  //JGTODO  TEST_EQUALITY(H.GetLevel(1)->Get< RCP<SmootherBase> >("PreSmoother", smooFactory)->GetType(),"Ifpack: Jacobi");
  //JGTODO  TEST_EQUALITY(H.GetLevel(1)->Get< RCP<SmootherBase> >("PostSmoother", smooFactory)->GetType(),"Ifpack: Jacobi");
#endif
    }
} //SetSmoothers

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver1)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");

  std::cout << "ICI" << std::endl;
  SmootherFactory SmooFactory(smoother);
  std::cout << "LA" << std::endl;

  RCP<Level> levelOne = rcp(new Level());

  Hierarchy H;
  H.SetLevel(levelOne);

  levelOne->Request("A");
  levelOne->Set("A",A);

  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  std::cout << "ICI2" << std::endl;
  H.SetCoarsestSolver(SmooFactory);
  std::cout << "LA2" << std::endl;

  RCP<SmootherBase>  preSmoo,postSmoo;
  preSmoo = levelOne->Get< RCP<SmootherBase> >("PreSmoother", &SmooFactory);
  TEST_INEQUALITY(preSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");
  postSmoo = levelOne->Get< RCP<SmootherBase> >("PostSmoother",&SmooFactory);

  TEST_INEQUALITY(postSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel");

#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver2)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  SmootherFactory SmooFactory(smoother);

  RCP<Level> levelOne = rcp(new Level());

  Hierarchy H;
  H.SetLevel(levelOne);
  levelOne->Request("A");
  levelOne->Set("A",A);

  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  H.SetCoarsestSolver(SmooFactory,MueLu::PRE);

  RCP<SmootherBase> preSmoo;
  preSmoo = levelOne->Get< RCP<SmootherBase> >("PreSmoother", &SmooFactory);

  TEST_INEQUALITY(preSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");

  TEST_EQUALITY(levelOne->IsAvailable("PostSmoother"), false);

#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver3)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  SmootherFactory SmooFactory(smoother);

  RCP<Level> levelOne = rcp(new Level());

  Hierarchy H;
  H.SetLevel(levelOne);
  levelOne->Request("A");
  levelOne->Set("A",A);

  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  H.SetCoarsestSolver(SmooFactory,MueLu::POST);

  TEST_EQUALITY(levelOne->IsAvailable("PreSmoother"), false);

  RCP<SmootherBase> postSmoo;
  postSmoo = levelOne->Get< RCP<SmootherBase> >("PostSmoother", &SmooFactory);

  TEST_INEQUALITY(postSmoo, Teuchos::null);
  //JG TODO  TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel");
#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_NoArgs)
{
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  H.SetLevel(levelOne);

  levelOne->Request("A");
  levelOne->Set("A",A);

  H.FullPopulate();
} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_AllArgs)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  H.SetLevel(levelOne);
  levelOne->Request("A");
  levelOne->Set("A",A);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory(PFact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));
  H.FullPopulate(PRFact,AcFact,SmooFact,0,2);
#endif
    }
} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,Iterate)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  H.SetLevel(Finest);

  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Set("NullSpace",nullSpace);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<CoalesceDropFactory> cdFact;
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<GenericPRFactory>   PRfact = rcp( new GenericPRFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );

#ifdef HAVE_MUELU_IFPACK
#ifdef HAVE_MUELU_AMESOS
  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel", 2);

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  int maxLevels = 5;
  status = H.FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(out,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  //Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  X->norm2(norms);
  X->scale(1/norms[0]);
  X->norm2(norms);
  out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  RHS->putScalar( (SC) 0.0);

  H.PrintResidualHistory(false);
  int iterations=10;
  H.Iterate(*RHS,iterations,*X);

  X->norm2(norms);
  out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(ios::fixed) <<
    std::setprecision(10) << norms[0] << std::endl;

  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
  TEST_EQUALITY(norms[0]<1e-10, true);

#endif
#endif
    }
} //Iterate

TEUCHOS_UNIT_TEST(Hierarchy,IterateWithImplicitRestriction)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.SetImplicitTranspose(true);
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  H.SetLevel(Finest);

  Finest->Request("A");
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);


  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);
  RCP<CoalesceDropFactory> cdFact;
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<GenericPRFactory>   PRfact = rcp( new GenericPRFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
  Acfact->SetImplicitTranspose(true);

#ifdef HAVE_MUELU_IFPACK
#ifdef HAVE_MUELU_AMESOS
  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel", 2);

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  int maxLevels = 5;
  status = H.FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(out,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  //Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  X->norm2(norms);
  X->scale(1/norms[0]);
  X->norm2(norms);
  out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  RHS->putScalar( (SC) 0.0);

  H.PrintResidualHistory(false);
  int iterations=10;
  H.Iterate(*RHS,iterations,*X);

  X->norm2(norms);
  out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(ios::fixed) <<
    std::setprecision(10) << norms[0] << std::endl;

  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
  TEST_EQUALITY(norms[0]<1e-10, true);

#endif
#endif
    }
} //Iterate

}//namespace MueLuTests

//Note from JG:
// For UnitTest,  TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Operator> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel");
// should be replaced by
// TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Operator> >("PreSmoother"), preSmoother);
// testing if preSmoother->GetType() == "Ifpack: Gauss-Seidel" should be a unit test of the class IfpackSmoother
