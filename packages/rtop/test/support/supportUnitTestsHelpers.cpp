
#include "supportUnitTestsHelpers.hpp"


int n = 4;

double errorTolSlack = 1e+1;

bool verbose = false;


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
  clp.setOption(
    "verbose", "quiet", &verbose,
    "Tests produce extra verbose output or not" );
}


} // namespace
