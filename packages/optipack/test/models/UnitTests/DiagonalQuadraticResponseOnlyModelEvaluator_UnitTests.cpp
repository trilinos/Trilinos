
#include "OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
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


template<class Scalar>
inline Scalar sqr(const Scalar &x) { return x*x; }


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DiagonalQuadraticResponseOnlyModelEvaluator,
  basic, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::derivativeGradient;
  using Thyra::create_DgDp_mv;
  using Thyra::eval_g_DgDp;
  using Thyra::get_mv;
  using Thyra::get_ele;
  using Thyra::norm_2;

  RCP<const Thyra::ModelEvaluator<Scalar> >
    model = OptiPack::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);
  
  TEST_ASSERT(!is_null(model));
  TEST_EQUALITY_CONST(model->Np(), 1);
  TEST_EQUALITY_CONST(model->Ng(), 1);

  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  RCP<const Thyra::VectorSpaceBase<Scalar> > g_space = model->get_g_space(0);

  RCP<Thyra::VectorBase<Scalar> > p_init = createMember(p_space);
  const Scalar val = as<Scalar>(2.0);
  out << "\nval = " << val << "\n";
  Thyra::V_S(p_init.ptr(), val);

  RCP<Thyra::VectorBase<Scalar> >
    g = createMember(g_space),
    g_grad = createMember(p_space);

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
    g.ptr(),derivativeGradient<Scalar>(g_grad) );

  out << "\ng =\n" << *g;
  out << "\ng_grad =\n" << *g_grad;

  const Ordinal globalDim = p_space->dim();
  out << "\nglobalDim = " << globalDim << "\n";
  
  TEST_FLOATING_EQUALITY( get_ele<Scalar>(*g, 0),
    as<Scalar>(0.5*globalDim)*val*val, as<ScalarMag>(g_tol/globalDim));
  
  TEST_FLOATING_EQUALITY(
    norm_2<Scalar>(*g_grad),
    ST::magnitude(ST::squareroot(as<Scalar>(globalDim)*val*val)),
    as<ScalarMag>(g_tol/globalDim));
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  DiagonalQuadraticResponseOnlyModelEvaluator, basic )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DiagonalQuadraticResponseOnlyModelEvaluator,
  offsets, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::VectorSpaceBase;
  using Thyra::VectorBase;
  using Thyra::derivativeGradient;
  using Thyra::create_DgDp_mv;
  using Thyra::eval_g_DgDp;
  using Thyra::get_mv;
  using Thyra::get_ele;
  using Thyra::norm_2;
  using Thyra::V_S;

  const RCP<OptiPack::DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
    model = OptiPack::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<const VectorSpaceBase<Scalar> > g_space = model->get_g_space(0);

  const Scalar p_soln_val = as<Scalar>(3.0);
  const RCP<VectorBase<Scalar> > p_soln = createMember(p_space);
  V_S(p_soln.ptr(), p_soln_val);
  model->setSolutionVector(p_soln);

  const Scalar g_offset = as<Scalar>(5.0);
  model->setScalarOffset(g_offset);

  const Scalar p_val = as<Scalar>(2.0);
  const RCP<VectorBase<Scalar> > p_init = createMember(p_space);
  V_S(p_init.ptr(), p_val);

  RCP<VectorBase<Scalar> >
    g = createMember(g_space),
    g_grad = createMember(p_space);

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
    g.ptr(),derivativeGradient<Scalar>(g_grad) );

  out << "\ng =\n" << *g;
  out << "\ng_grad =\n" << *g_grad;

  const Ordinal globalDim = p_space->dim();
  out << "\nglobalDim = " << globalDim << "\n";
  
  TEST_FLOATING_EQUALITY(
    get_ele<Scalar>(*g, 0),
    as<Scalar>(0.5 * sqr(p_val - p_soln_val) * globalDim + g_offset),
    as<ScalarMag>(g_tol/globalDim)
    );
  
  TEST_FLOATING_EQUALITY(
    sum(*g_grad),
    as<Scalar>( (p_val - p_soln_val) * globalDim ),
    as<ScalarMag>(g_tol/globalDim)
    );
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  DiagonalQuadraticResponseOnlyModelEvaluator, offsets )


} // namespace


