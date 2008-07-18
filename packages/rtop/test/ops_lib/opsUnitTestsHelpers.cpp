
#include "opsUnitTestsHelpers.hpp"


int n = 4;

double errorTolSlack = 1e+1;


namespace {


TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp =
    Teuchos::UnitTestRepository::getCLP();
  clp.setOption(
    "n", &n, "Number of elements in the local vectors" );
  clp.setOption(
    "error-tol-slack", &errorTolSlack,
    "Slack off of machine epsilon used to check test results" );
}


} // namespace
