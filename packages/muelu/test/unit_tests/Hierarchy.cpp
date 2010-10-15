#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Hierarchy.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Hierarchy,Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  typedef Tpetra::Map<LO,GO,Node> Map;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> Operator;
  typedef Tpetra::Vector<Scalar,LO,GO,Node>    Vector;
  typedef MueLu::Level<Scalar,LO,GO,Node>    Level;

  typedef MueLu::Hierarchy<Scalar,LO,GO,Node>    Hierarchy;

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  Level levelOne;
  levelOne.SetLevelID(1);
  Level levelTwo;
  levelTwo.SetLevelID(2);

  Hierarchy H;

  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);

  std::cout << H << std::endl;

  /* TODO
  Test set/get of R & P matrices.
  */
  /*
  RCP<Operator> x = rcp(new Tpetra::Vector<Scalar,LO,GO,Node>(map,nx) );
  RCP<Operator> y = rcp(new Tpetra::Vector<Scalar,LO,GO,Node>(map,nx) );
  x->putScalar(1);
  */


}

}//namespace <anonymous>

