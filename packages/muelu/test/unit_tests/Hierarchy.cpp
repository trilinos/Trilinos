#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Hierarchy,Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

/*
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef Cthulhu::Operator<Scalar,LO,GO,Node> Operator;
  typedef Cthulhu::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node,LMO>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node,LMO>    Hierarchy;
*/

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<Level> levelTwo = rcp(new Level() );
  levelTwo->SetLevelID(2);

  Hierarchy H;

  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);

  std::cout << H << std::endl;

  /* TODO
  Test set/get of R & P matrices.
  */
  /*
  RCP<Vector> x = rcp(new Cthulhu::TpetraVector<Scalar,LO,GO,Node>(map,nx) );
  RCP<Vector> y = rcp(new Cthulhu::TpetraVector<Scalar,LO,GO,Node>(map,nx) );
  x->putScalar(1);
  */

} //TEST0

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy1)
{

/*
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef Cthulhu::Operator<Scalar,LO,GO,Node> Operator;
  typedef Cthulhu::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node,LMO>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node,LMO>    Hierarchy;
*/

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  try {
    out << "Intentionally providing no prolongator factory to FillHierarchy .... ";

    RCP< SaPFactory >    PFact = Teuchos::null;
    //RCP<OperatorFactory<Scalar,LO,GO,Node> >  opFact = PFact;
    //H.FillHierarchy(opFact);
    H.FillHierarchy(PFact);
    //H.FillHierarchy(Teuchos::null);
  }
  catch(...) {
    out << "Caught the error" << std::endl;
  }

}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy2)
{

/*
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef Cthulhu::Operator<Scalar,LO,GO,Node> Operator;
  typedef Cthulhu::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node,LMO>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node,LMO>    Hierarchy;
*/

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory<Scalar,LO,GO,Node,LMO> >    PFact = rcp(new SaPFactory<Scalar,LO,GO,Node,LMO>());

  out << "Providing just prolongator factory to FillHierarchy." << std::endl;
  H.FillHierarchy(PFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy3)
{

/*
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef Cthulhu::Operator<Scalar,LO,GO,Node> Operator;
  typedef Cthulhu::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node,LMO>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node,LMO>    Hierarchy;
*/

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory<Scalar,LO,GO,Node,LMO> >    PFact = rcp(new SaPFactory<Scalar,LO,GO,Node,LMO>());
  RCP<TransPFactory<Scalar,LO,GO,Node,LMO> > RFact = rcp(new TransPFactory<Scalar,LO,GO,Node,LMO>());
  RCP<RAPFactory<Scalar,LO,GO,Node,LMO> >    AcFact= rcp(new RAPFactory<Scalar,LO,GO,Node,LMO>());

  out << "Providing all three factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PFact,RFact,AcFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{

/*
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef Cthulhu::Operator<Scalar,LO,GO,Node> Operator;
  typedef Cthulhu::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node,LMO>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node,LMO>    Hierarchy;
*/

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);
  H.SetSmoothers();
}

}//namespace <anonymous>

