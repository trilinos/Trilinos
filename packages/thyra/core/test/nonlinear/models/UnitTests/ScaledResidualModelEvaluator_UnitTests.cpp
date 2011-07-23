#include "Thyra_ScaledResidualModelEvaluator.hpp"
#include "Thyra_Simple2DModelEvaluator.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {

using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Thyra::createMember;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ScalarResidualModelEvaluator,
  basic, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::derivativeGradient;

  RCP<Thyra::ModelEvaluator<Scalar> > model = 
    Thyra::simple2DModelEvaluator<Scalar>();
  
  RCP<Thyra::ScaledResidualModelEvaluator<Scalar> > scaled_model = 
    rcp(new Thyra::ScaledResidualModelEvaluator<Scalar>(model));
  
  TEST_ASSERT(!is_null(model));

  Thyra::ModelEvaluatorBase::InArgs<Scalar> in_args = model->getNominalValues();

  RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(model->get_x_space());
  Thyra::V_S(x.ptr(), Teuchos::as<Scalar>(2.0));
  in_args.set_x(x);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> out_args = model->createOutArgs();

  RCP<const Thyra::VectorSpaceBase<Scalar> > f_space = model->get_f_space();

  RCP<Thyra::VectorBase<Scalar> > f = createMember(f_space);
  RCP<Thyra::VectorBase<Scalar> > f_scaled = createMember(f_space);
  
  RCP<Thyra::LinearOpBase<Scalar> > W = model->create_W_op() ;
  RCP<Thyra::LinearOpBase<Scalar> > W_scaled = model->create_W_op() ;

  RCP<Thyra::VectorBase<Scalar> > scaling_diagonal = createMember(f_space);
  const Scalar val = as<Scalar>(2.0);
  Thyra::V_S(scaling_diagonal.ptr(), val);
  
  scaled_model->setScalingDiagonal(scaling_diagonal);

  out_args.set_f(f);
  out_args.set_W_op(W);

  model->evalModel(in_args, out_args);

  out_args.set_f(f_scaled);
  out_args.set_W_op(W_scaled);

  scaled_model->evalModel(in_args, out_args);

  ScalarMag tol = Teuchos::as<ScalarMag>(10.0) * ST::eps();

  TEST_FLOATING_EQUALITY(2.0 * Thyra::get_ele(*f,0), Thyra::get_ele(*f_scaled,0),
			 tol);
  TEST_FLOATING_EQUALITY(2.0 * Thyra::get_ele(*f,1), Thyra::get_ele(*f_scaled,1),
			 tol);

  RCP<Thyra::MultiVectorBase<Scalar> > M = 
    Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W);
  TEUCHOS_ASSERT(Teuchos::nonnull(M));
  Thyra::DetachedMultiVectorView<Scalar> M_dv(*M);

  RCP<Thyra::MultiVectorBase<Scalar> > M_scaled = 
    Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_scaled);
  TEUCHOS_ASSERT(Teuchos::nonnull(M_scaled));
  Thyra::DetachedMultiVectorView<Scalar> M_dv_scaled(*M_scaled);


  TEST_FLOATING_EQUALITY(2.0 * M_dv(0,0), M_dv_scaled(0,0), tol);
  TEST_FLOATING_EQUALITY(2.0 * M_dv(0,1), M_dv_scaled(0,1), tol);
  TEST_FLOATING_EQUALITY(2.0 * M_dv(1,0), M_dv_scaled(1,0), tol);
  TEST_FLOATING_EQUALITY(2.0 * M_dv(1,1), M_dv_scaled(1,1), tol);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  ScalarResidualModelEvaluator, basic )

} // namespace


