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

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;

TEUCHOS_UNIT_TEST(Hierarchy,Constructor)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEUCHOS_TEST_INEQUALITY(H, Teuchos::null, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy,SetAndGetLevel)
{


  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> level = rcp(new Level() );
  H.SetLevel(level);
  RCP<Level> dupLevel = H.GetLevel(0);
  TEUCHOS_TEST_EQUALITY(level, dupLevel, out, success);

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
  TEUCHOS_TEST_EQUALITY(H.GetNumberOfLevels(), 3, out, success);

}//NumberOfLevels

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_NoFactoriesGiven)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

  out << "Providing no factories to FillHierarchy." << std::endl;
  H.FillHierarchy();

  bool ToF = H.PrintResidualHistory();
  H.PrintResidualHistory(!ToF);
  TEUCHOS_TEST_INEQUALITY(H.PrintResidualHistory(), ToF, out, success);
} //FillHierarchy_NoFactoriesGiven

#ifdef FAILING_ONE_PROC_TERMINATING_BADLY_DONT_KNOW_WHY
TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PRFactoryOnly)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just PR factory to FillHierarchy." << std::endl;
  Teuchos::ParameterList status;
  status = H.FillHierarchy(PRFact);
  TEUCHOS_TEST_EQUALITY(status.get("fine nnz",(Xpetra::global_size_t)-1), 295, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("total nnz",(Xpetra::global_size_t)-1), 422, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("start level",-1), 0, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("end level",-1), 1, out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(status.get("operator complexity",(SC)-1.0),1.43051,1e-5,out,success);

} //FillHierarchy_PRFactoryOnly
#endif

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_BothFactories)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<TransPFactory>  RFact = rcp(new TransPFactory());
  GenericPRFactory PRFact(PFact,RFact);
  RAPFactory    AcFact;

  out << "Providing both factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact,AcFact);
} //FillHierarchy_BothFactories

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {


  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

//   H.SetSmoothers();

//   TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->template Get< Teuchos::RCP<SmootherPrototype> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel", out, success);
//   TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->template Get< Teuchos::RCP<SmootherPrototype> >("PostSmoother")->GetType(),"Ifpack: Gauss-Seidel", out, success);

#ifdef HAVE_MUELU_IFPACK
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Jacobi");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto) );
  H.SetSmoothers(*smooFactory);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->Get< Teuchos::RCP<SmootherPrototype> >("PreSmoother")->GetType(),"Ifpack: Jacobi", out, success);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->Get< Teuchos::RCP<SmootherPrototype> >("PostSmoother")->GetType(),"Ifpack: Jacobi", out, success);
#endif
    }
} //SetSmoothers

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

#ifdef HAVE_MUELU_IFPACK
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype> smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  SmootherFactory SmooFactory(smoother);

  H.SetCoarsestSolver(SmooFactory);

  RCP<SmootherPrototype>  preSmoo,postSmoo;

  preSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PreSmoother");
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PostSmoother");
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::PRE);
  preSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PreSmoother");
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PostSmoother");
  TEUCHOS_TEST_EQUALITY(postSmoo, Teuchos::null, out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::POST);
  preSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PreSmoother");
  TEUCHOS_TEST_EQUALITY(preSmoo, Teuchos::null, out, success);
  postSmoo = levelOne->Get< Teuchos::RCP<SmootherPrototype> >("PostSmoother");
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_NoArgs)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

  H.FullPopulate();
    }
} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_AllArgs)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());
  levelOne->Set("A",A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory(PFact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory());

#ifdef HAVE_MUELU_IFPACK
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));
  H.FullPopulate(PRFact,AcFact,SmooFact,0,2);
#endif
    }
} //FullPopulate

#include "MueLu_AggregationOptions.hpp"

TEUCHOS_UNIT_TEST(Hierarchy,Iterate)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Set("NullSpace",nullSpace);
  H.SetLevel(Finest);

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
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 2);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

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
  TEUCHOS_TEST_EQUALITY(norms[0]<1e-10, true, out, success);

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
  RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
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

  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Set("NullSpace",nullSpace);
  H.SetLevel(Finest);

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
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 2);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

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
  TEUCHOS_TEST_EQUALITY(norms[0]<1e-10, true, out, success);

#endif
#endif
    }
} //Iterate

}//namespace <anonymous>

//Note from JG:
// For UnitTest,  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->Get< Teuchos::RCP<Operator> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel", out, success);
// should be replaced by
// TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->Get< Teuchos::RCP<Operator> >("PreSmoother"), preSmoother, out, success);
// testing if preSmoother->GetType() == "Ifpack: Gauss-Seidel" should be a unit test of the class IfpackSmoother
