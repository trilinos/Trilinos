#include "Teuchos_UnitTestHarness.hpp"
//#include "test_helpers.hpp"

#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor)
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

  typedef MueLu::SmootherFactory<Scalar,LO,GO,Node,LMO> SmootherFactory;
*/

  using namespace Teuchos;
  //using namespace MueLu;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SmootherPrototype>  smoother = Teuchos::null;
  RCP<SmootherFactory> SmooFactory = rcp(new SmootherFactory(smoother) );
  TEUCHOS_TEST_EQUALITY(SmooFactory != Teuchos::null, true, out, success);

}

}//namespace <anonymous>

