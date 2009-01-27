
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DiagonalQuadraticResponseOnlyModelEvaluator, basic, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::create_DgDp_mv;
  using Thyra::eval_g_DgDp;
  using Thyra::get_ele;
  using Thyra::norm_2;

  ECHO(RCP<const Thyra::ModelEvaluator<Scalar> >
    model = OptiPack::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim));
  
  TEST_ASSERT(!is_null(model));
  TEST_EQUALITY_CONST(model->Np(), 1);
  TEST_EQUALITY_CONST(model->Ng(), 1);

  ECHO(RCP<const Thyra::VectorSpaceBase<Scalar> > p_space = model->get_p_space(0));
  ECHO(RCP<const Thyra::VectorSpaceBase<Scalar> > g_space = model->get_g_space(0));

  ECHO(RCP<Thyra::VectorBase<Scalar> > p_init = createMember(p_space));
  ECHO(const Scalar val = as<Scalar>(2.0));
  out << "\nval = " << val << "\n";
  ECHO(Thyra::V_S(p_init.ptr(), val));

  ECHO(RCP<Thyra::VectorBase<Scalar> > g = createMember(g_space));
  ECHO(MEB::Derivative<Scalar> DgDp =
    create_DgDp_mv<Scalar>(*model, 0, 0, MEB::DERIV_TRANS_MV_BY_ROW));

  ECHO(eval_g_DgDp<Scalar>(*model, 0, *p_init, 0, g.ptr(), DgDp));

  out << "\ng =\n" << *g;
  out << "\nDgDp =\n" << *DgDp.getDerivativeMultiVector().getMultiVector();

  ECHO(const Ordinal globalDim = p_space->dim());
  out << "\nglobalDim = " << globalDim << "\n";
  
  TEST_FLOATING_EQUALITY( get_ele<Scalar>(*g, 0),
    as<Scalar>(0.5*globalDim)*val*val, as<ScalarMag>(g_tol/globalDim));
  
  TEST_FLOATING_EQUALITY(
    norm_2<Scalar>(*DgDp.getDerivativeMultiVector().getMultiVector()->col(0)),
    ST::squareroot(as<Scalar>(globalDim)*val*val), as<ScalarMag>(g_tol/globalDim));
    
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DiagonalQuadraticResponseOnlyModelEvaluator, basic )


} // namespace


