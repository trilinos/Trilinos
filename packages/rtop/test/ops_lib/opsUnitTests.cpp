#include "RTOpPack_ROpDotProd.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_as.hpp"


namespace {


// Size of the vectors
int n = 4;
// Cushion off of machine eps
double errorTolSlack = 1e-2;


class UnitTestSetup {
public:
  UnitTestSetup()
    {
      Teuchos::CommandLineProcessor &clp =
        Teuchos::UnitTestRepository::getCLP();
      clp.setOption(
        "n", &n, "Number of elements in the local vectors" );
      clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
    }
} unitTestSetup;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DotProd, basic, T )
{
  out << "\nToDo: Implement!\n";
  success = false;
}

//
// Test instantations
//


#define UNIT_TEST_GROUP(T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DotProd, basic, T )


UNIT_TEST_GROUP( float )
UNIT_TEST_GROUP( double )

} // namespace
