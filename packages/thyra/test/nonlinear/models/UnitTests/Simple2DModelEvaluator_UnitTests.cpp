
#include "Thyra_Simple2DModelEvaluator.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
typedef Thyra::ModelEvaluatorBase MEB;
using Thyra::Simple2DModelEvaluator;
using Thyra::simple2DModelEvaluator;


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleModelEvaluator, construct, Scalar )
{
  RCP<Simple2DModelEvaluator<Scalar> > model = simple2DModelEvaluator<Scalar>();
  TEST_ASSERT(model != null);
  TEST_EQUALITY(model->Np(), 0);
  TEST_EQUALITY(model->Ng(), 0);
  TEST_ASSERT(model->get_x_space() != null);
  TEST_EQUALITY(model->get_x_space()->dim(), 2);
  TEST_ASSERT(model->get_f_space() != null);
  TEST_EQUALITY(model->get_f_space()->dim(), 2);
  // ToDo: Test getNominalValues()
  TEST_ASSERT(model->create_W_op() != null);
  TEST_ASSERT(model->get_W_factory() != null);
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_EQUALITY(inArgs.Np(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  SimpleModelEvaluator, construct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleModelEvaluator, eval, Scalar )
{
  RCP<Simple2DModelEvaluator<Scalar> > model = simple2DModelEvaluator<Scalar>();
  // ToDo: Finish this!
  //TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  SimpleModelEvaluator, eval )


} // namespace
