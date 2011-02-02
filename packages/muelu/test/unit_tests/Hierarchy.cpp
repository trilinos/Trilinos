#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PRFactory.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

TEUCHOS_UNIT_TEST(Hierarchy,Constructor)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEUCHOS_TEST_INEQUALITY(H, Teuchos::null, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy,SetAndGetLevel)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> level = rcp(new Level() );
  H.SetLevel(level);
  RCP<Level> dupLevel = H.GetLevel(0);
  TEUCHOS_TEST_EQUALITY(level, dupLevel, out, success);

}//SetAndGetLevel

TEUCHOS_UNIT_TEST(Hierarchy,NumberOfLevels)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

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
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  out << "Providing no factories to FillHierarchy." << std::endl;
  H.FillHierarchy();

  bool ToF = H.PrintResidualHistory();
  H.PrintResidualHistory(!ToF);
  TEUCHOS_TEST_INEQUALITY(H.PrintResidualHistory(), ToF, out, success);
} //FillHierarchy_NoFactoriesGiven

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PRFactoryOnly)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just PR factory to FillHierarchy." << std::endl;
  Teuchos::ParameterList status;
  status = H.FillHierarchy(PRFact);
  TEUCHOS_TEST_EQUALITY(status.get("fine nnz",(Cthulhu::global_size_t)-1), 295, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("total nnz",(Cthulhu::global_size_t)-1), 392, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("start level",-1), 0, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("end level",-1), 1, out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(status.get("operator complexity",(SC)-1.0),1.32881,1e-5,out,success);

} //FillHierarchy_PRFactoryOnly

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_BothFactories)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);
  RAPFactory    AcFact;

  out << "Providing both factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact,AcFact);
} //FillHierarchy_BothFactories

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);
  H.SetSmoothers();

  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPreSmoother()->GetType(),"Ifpack: Gauss-Seidel", out, success);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPostSmoother()->GetType(),"Ifpack: Gauss-Seidel", out, success);


  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Jacobi");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto) );
  H.SetSmoothers(smooFactory);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPreSmoother()->GetType(),"Ifpack: Jacobi", out, success);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPostSmoother()->GetType(),"Ifpack: Jacobi", out, success);

} //SetSmoothers

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  SmootherFactory SmooFactory(smoother);

  H.SetCoarsestSolver(SmooFactory);

  RCP<SmootherPrototype>  preSmoo,postSmoo;

  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::PRE);
  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_EQUALITY(postSmoo, Teuchos::null, out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::POST);
  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_EQUALITY(preSmoo, Teuchos::null, out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);

} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_NoArgs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  H.FullPopulate();

} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_AllArgs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory(PFact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory());

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));

  H.FullPopulate(PRFact,AcFact,SmooFact,0,2);

} //FullPopulate

}//namespace <anonymous>
