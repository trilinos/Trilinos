
#include "OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Thyra::createMember;


Teuchos_Ordinal g_localDim = 4;

double g_tol = Teuchos::ScalarTraits<double>::eps();


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Number of local vector elements on each process" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol, "Floating point tolerance" );
}


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, basic, Scalar )
{

  //TEST_FOR_EXCEPT(true);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( NonlinearCG, basic )


} // namespace


