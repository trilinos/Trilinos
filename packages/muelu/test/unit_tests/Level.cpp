#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Level.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Level, Test0)
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

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  Teuchos::ParameterList list;
  list.set("nx",10);
  LO nx = list.get("nx",10);
  RCP<const Map> map = MueLu_UnitTest::create_tpetra_map<LO,GO,Node>(nx);
  RCP<Operator> A = CreateCrsMatrix<Scalar,LO,GO,Node>("Laplace1D",map,list);

  out << "Testing default ctor" << std::endl;
  Level firstLevel;
  out << "Testing set methods" << std::endl;
  firstLevel.SetLevelID(1);
  firstLevel.SetA(A);

  /* TODO
  Test set/get of R & P matrices.
  */
  /*
  RCP<Operator> x = rcp(new Tpetra::Vector<Scalar,LO,GO,Node>(map,nx) );
  RCP<Operator> y = rcp(new Tpetra::Vector<Scalar,LO,GO,Node>(map,nx) );
  x->putScalar(1);
  */


  //out << firstLevel << std::endl;
  out << "Testing copy ctor" << std::endl;
  Level secondLevel(firstLevel);
  //out << secondLevel << std::endl;
}

}//namespace <anonymous>

