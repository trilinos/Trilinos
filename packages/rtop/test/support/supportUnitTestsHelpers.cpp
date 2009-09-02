
#include "supportUnitTestsHelpers.hpp"


int TestingSupportHelpers::n = 4;

double TestingSupportHelpers::errorTolSlack = 1e+1;

bool TestingSupportHelpers::verbose = false;


namespace {


using TestingSupportHelpers::n;
using TestingSupportHelpers::errorTolSlack;
using TestingSupportHelpers::verbose;


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
