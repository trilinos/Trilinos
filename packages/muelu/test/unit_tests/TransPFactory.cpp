#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TransPFactory.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(TransPFactory, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef double ScalarType;
  typedef double Scalar;
  typedef int    LO;
  typedef int    GO;
  typedef Kokkos::DefaultNode::DefaultNodeType Node;
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO;

  typedef MueLu::TransPFactory<Scalar,LO,GO,Node,LMO>    TransPFactory;

  using namespace Teuchos;
  using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<TransPFactory> transPFactory = rcp(new TransPFactory);
  TEUCHOS_TEST_EQUALITY(transPFactory != Teuchos::null, true, out, success);

  out << *transPFactory << std::endl;

}

}//namespace <anonymous>

