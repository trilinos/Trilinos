
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


const Teuchos_Ordinal g_localDim = 4; // ToDo: Make variable!


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DiagonalQuadraticResponseOnlyModelEvaluator, basic, Scalar )
{
  RCP<Thyra::ModelEvaluator<Scalar> >
    model = OptiPack::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);
  
  TEST_ASSERT(!is_null(model));
  TEST_EQUALITY_CONST(model->Np(), 1);
  TEST_EQUALITY_CONST(model->Ng(), 1);
  ECHO(RCP<const Thyra::VectorSpaceBase<Scalar> > p_space = model->get_p_space(0));
  ECHO(RCP<const Thyra::VectorSpaceBase<Scalar> > g_space = model->get_g_space(0));
  ECHO(RCP<Thyra::VectorBase<Scalar> > p_init = createMember(p_space));
  ECHO(Thyra::V_S(p_init.ptr(), as<Scalar>(1.0)));
  ECHO(RCP<Thyra::VectorBase<Scalar> > g = createMember(g_space));
  ECHO(Thyra::eval_g<Scalar>(*model, 0, *p_init, 0, &*g));
  out << "\ng =\n" << *g;
  
}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DiagonalQuadraticResponseOnlyModelEvaluator, basic )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DiagonalQuadraticResponseOnlyModelEvaluator, basic, double )


} // namespace


