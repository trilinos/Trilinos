
#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_TestingTools.hpp"

#include "OptiPack_NonlinearCG.hpp"
#include "OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "GlobiPack_ArmijoPolyInterpLineSearch.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Thyra::createMember;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using GlobiPack::ArmijoPolyInterpLineSearch;
using GlobiPack::armijoQuadraticLineSearch;
using OptiPack::NonlinearCG;
using OptiPack::nonlinearCG;
using OptiPack::DiagonalQuadraticResponseOnlyModelEvaluator;
using OptiPack::diagonalQuadraticResponseOnlyModelEvaluator;
namespace NCGU = OptiPack::NonlinearCGUtils;


Teuchos_Ordinal g_localDim = 4;

double g_tol_scale = 100.0;

double g_p_init = 1.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Number of local vector elements on each process" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol-scale", &g_tol_scale, "Floating point tolerance" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "p-init", &g_p_init, "Initial guess for unknowns p" );
}


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_basic, Scalar )
{

  using Teuchos::optInArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<VectorBase<Scalar> > ps = createMember(p_space);
  const Scalar ps_val = 2.0;
  V_S(ps.ptr(), ps_val);
  model->setSolutionVector(ps);
  const ScalarMag g_offset = as<Scalar>(5.0);
  model->setScalarOffset(g_offset);

  // Set up a quadratic interploation line search that will do just one
  // iteration and should exactly minimize a quadratic function.
  const RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();
  const RCP<ParameterList> lsPL = parameterList();
  lsPL->set("Min Backtrack Fraction", 0.0);
  lsPL->set("Max Backtrack Fraction", 1e+50);
  lsPL->set("Min Number of Iterations", 1);
  lsPL->set("Max Number of Iterations", 2);
  linesearch->setParameterList(lsPL);

  const RCP<NonlinearCG<Scalar> > cgSolver =
    nonlinearCG<Scalar>(model, 0, 0, linesearch);

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  const ScalarMag tol = as<Scalar>(g_tol_scale) * ST::eps();
  const ScalarMag alpha_init = 10.0;
  ScalarMag g_opt = -1.0;
    
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(alpha_init) );
 
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_FLOATING_EQUALITY(g_opt, g_offset, tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *ps,
    "tol", tol,
    "2*tol", as<ScalarMag>(2.0)*tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, FR_basic )


} // namespace


