
#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_TestingTools.hpp"

#include "OptiPack_NonlinearCG.hpp"
#include "OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
#include "GlobiPack_ArmijoPolyInterpLineSearch.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ArrayView;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::ScalarTraits;
using Teuchos::FancyOStream;
using Thyra::createMember;
using RTOpPack::ReductTarget;
using RTOpPack::ConstSubVectorView;
using RTOpPack::SubVectorView;
typedef Thyra::Ordinal Ordinal;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::applyOp;
using GlobiPack::ArmijoPolyInterpLineSearch;
using GlobiPack::armijoQuadraticLineSearch;
using OptiPack::NonlinearCG;
using OptiPack::nonlinearCG;
using OptiPack::DiagonalQuadraticResponseOnlyModelEvaluator;
using OptiPack::diagonalQuadraticResponseOnlyModelEvaluator;
namespace NCGU = OptiPack::NonlinearCGUtils;


Teuchos_Ordinal g_globalDim = 16;

double g_solve_tol_scale = 10.0;

double g_error_tol_scale = 1000.0;

double g_nonlin_term_factor = 1e-2;

double g_nonlin_solve_tol = 1e-4;

double g_nonlin_error_tol = 1e-3;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "global-dim", &g_globalDim,
    "Number of global vector elements over all processes" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "solve-tol-scale", &g_solve_tol_scale,
    "Floating point tolerance for nonlinear CG solve for linear CG tests" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "error-tol-scale", &g_error_tol_scale,
    "Floating point tolerance for error checks for linear CG tests" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "nonlin-term-factor", &g_nonlin_term_factor,
    "Scale factor for cubic term in objective" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "nonlin-solve-tol", &g_nonlin_solve_tol,
    "Floating point tolerance for general nonlinear CG solve" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "nonlin-error-tol", &g_nonlin_error_tol,
    "Floating point tolerance for error checks for general nonlinear CG solve" );

}


template<class Scalar>
const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
createModel(
  const int globalDim,
  const typename ScalarTraits<Scalar>::magnitudeType & g_offset
  )
{

  const RCP<const Teuchos::Comm<Thyra::Ordinal> > comm =
    Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

  const int numProcs = comm->getSize();
  TEST_FOR_EXCEPT_MSG( numProcs > globalDim,
    "Error, the number of processors can not be greater than the global"
    " dimension of the vectors!." );
  const int localDim = globalDim / numProcs;
  const int localDimRemainder = globalDim % numProcs;
  TEST_FOR_EXCEPT_MSG( localDimRemainder != 0,
    "Error, the number of processors must divide into the global number"
    " of elements exactly for now!." );

  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(localDim);
  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<VectorBase<Scalar> > ps = createMember(p_space);
  const Scalar ps_val = 2.0;
  V_S(ps.ptr(), ps_val);
  model->setSolutionVector(ps);
  model->setScalarOffset(g_offset);

  return model;

}


template<class Scalar>
const RCP<NonlinearCG<Scalar> >
createNonlinearCGSolver(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const RCP<FancyOStream> &out
  )
{
 
  // Set up a quadratic interploation line search that will do just one
  // iteration and should exactly minimize a quadratic function.
  const RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();
  const RCP<ParameterList> lsPL = parameterList();
  lsPL->set("Armijo Slope Fraction", 0.0);
  lsPL->set("Min Backtrack Fraction", 0.0);
  lsPL->set("Max Backtrack Fraction", 1e+50);
  lsPL->set("Min Num Iterations", 1);
  lsPL->set("Max Num Iterations", 2);
  linesearch->setParameterList(lsPL);

  const RCP<NonlinearCG<Scalar> > cgSolver =
    nonlinearCG<Scalar>(model, 0, 0, linesearch);

  const RCP<ParameterList> pl = parameterList();
  //pl->set("AND Convergence Tests", true);
  cgSolver->setParameterList(pl);

  cgSolver->setOStream(out);

  return cgSolver;

}


//
// RTOp to Assign elements z[i] = i + 1, i = 0...n-1
//

template<class Scalar>
class TOpAssignValToGlobalIndex : public RTOpPack::RTOpT<Scalar> {
public:
  TOpAssignValToGlobalIndex() {}
protected:
  bool coord_invariant_impl() const
    {
      return true;
    }
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;
      validate_apply_op( *this, 0, 1, false, sub_vecs, targ_sub_vecs, reduct_obj );
      const SubVectorView<Scalar> &z = targ_sub_vecs[0];
      const Ordinal z_global_offset = z.globalOffset();
      const Ordinal z_sub_dim = z.subDim();
      iter_t z_val = z.values().begin();
      const ptrdiff_t z_val_s = z.stride();

      for ( int i = 0; i < z_sub_dim; ++i, z_val += z_val_s ) {
        *z_val = as<Scalar>(z_global_offset + i + 1);
      }
    }
};



//
// Unit tests
//


//
// Check that internal default parameters are set correctly
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, defaultParams, Scalar )
{
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  namespace NCGU = OptiPack::NonlinearCGUtils;
  const RCP<NonlinearCG<Scalar> > cgSolver = nonlinearCG<Scalar>();
  TEST_EQUALITY(cgSolver->get_alpha_init(), as<ScalarMag>(NCGU::alpha_init_default));
  TEST_EQUALITY(cgSolver->get_alpha_reinit(), NCGU::alpha_reinit_default);
  TEST_EQUALITY(cgSolver->get_minIters(), NCGU::minIters_default);
  TEST_EQUALITY(cgSolver->get_maxIters(), NCGU::maxIters_default);
  TEST_EQUALITY(cgSolver->get_g_reduct_tol(), as<ScalarMag>(NCGU::g_reduct_tol_default));
  TEST_EQUALITY(cgSolver->get_g_grad_tol(), as<ScalarMag>(NCGU::g_grad_tol_default));
  TEST_EQUALITY(cgSolver->get_g_mag(), as<ScalarMag>(NCGU::g_mag_default));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, defaultParams )


//
// Check that internal default parameters are set correctly when gotten off of
// an empty parameter list
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, parseParamsDefaultParams, Scalar )
{
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  namespace NCGU = OptiPack::NonlinearCGUtils;
  const RCP<NonlinearCG<Scalar> > cgSolver = nonlinearCG<Scalar>();
  const RCP<ParameterList> pl = parameterList();
  cgSolver->setParameterList(pl);
  TEST_EQUALITY(cgSolver->get_alpha_init(), as<ScalarMag>(NCGU::alpha_init_default));
  TEST_EQUALITY(cgSolver->get_alpha_reinit(), NCGU::alpha_reinit_default);
  TEST_EQUALITY(cgSolver->get_minIters(), NCGU::minIters_default);
  TEST_EQUALITY(cgSolver->get_maxIters(), NCGU::maxIters_default);
  TEST_EQUALITY(cgSolver->get_g_reduct_tol(), as<ScalarMag>(NCGU::g_reduct_tol_default));
  TEST_EQUALITY(cgSolver->get_g_grad_tol(), as<ScalarMag>(NCGU::g_grad_tol_default));
  TEST_EQUALITY(cgSolver->get_g_mag(), as<ScalarMag>(NCGU::g_mag_default));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, parseParamsDefaultParams )


//
// Check that parameter list is parsed correctly
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, parseParams, Scalar )
{
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  namespace NCGU = OptiPack::NonlinearCGUtils;
  const RCP<NonlinearCG<Scalar> > cgSolver = nonlinearCG<Scalar>();
  const double alpha_init = 0.9;
  const bool alpha_reinit = true;
  const int minIters = 92;
  const int maxIters = 99;
  const double g_reduct_tol = 2.5;
  const double g_grad_tol = 2.8;
  const double g_mag = 3.1;
  TEST_INEQUALITY( alpha_reinit, NCGU::alpha_reinit_default ); // Make sure different
  const RCP<ParameterList> pl = parameterList();
  pl->set("Initial Linesearch Step Length", alpha_init);
  pl->set("Reinitlaize Linesearch Step Length", alpha_reinit);
  pl->set("Min Num Iterations", minIters);
  pl->set("Max Num Iterations", maxIters);
  pl->set("Objective Reduction Tol", g_reduct_tol);
  pl->set("Objective Gradient Tol", g_grad_tol);
  pl->set("Objective Magnitude", g_mag);
  cgSolver->setParameterList(pl);
  const ScalarMag tol = SMT::eps();
  TEST_FLOATING_EQUALITY(cgSolver->get_alpha_init(), as<ScalarMag>(alpha_init), tol);
  TEST_EQUALITY(cgSolver->get_alpha_reinit(), alpha_reinit);
  TEST_EQUALITY(cgSolver->get_minIters(), minIters);
  TEST_EQUALITY(cgSolver->get_maxIters(), maxIters);
  TEST_FLOATING_EQUALITY(cgSolver->get_g_reduct_tol(), as<ScalarMag>(g_reduct_tol), tol);
  TEST_FLOATING_EQUALITY(cgSolver->get_g_grad_tol(), as<ScalarMag>(g_grad_tol), tol);
  TEST_FLOATING_EQUALITY(cgSolver->get_g_mag(), as<ScalarMag>(g_mag), tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, parseParams )


//
// Print valid parameters
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, printValidParams, Scalar )
{
  const RCP<NonlinearCG<Scalar> > cgSolver = nonlinearCG<Scalar>();
  const RCP<const ParameterList> validPL = cgSolver->getValidParameters();
  typedef Teuchos::ParameterList::PrintOptions PO;
  out << "\nvalidPL:\n";
  validPL->print(out, PO().indent(2).showTypes(true).showFlags(2).showDoc(2));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, printValidParams )


//
// Test basic convergence in one iteration for one eignvalue
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_oneEigenVal, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);
  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const int dim = p_space->dim();

  const RCP<NonlinearCG<Scalar> > cgSolver =
    createNonlinearCGSolver<Scalar>(model, rcpFromRef(out));

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  ScalarMag g_opt = -1.0;
  const ScalarMag tol = as<Scalar>(g_solve_tol_scale * dim) * ST::eps();
  const ScalarMag alpha_init = 10.0;
  int numIters = -1;
    
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(tol), optInArg(alpha_init), outArg(numIters) );

  out << "\n";
 
  const ScalarMag err_tol = as<Scalar>(g_error_tol_scale * dim) * ST::eps();
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_EQUALITY( numIters, 1);
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, FR_oneEigenVal )


//
// Test convergence for partially unique eigen values
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_partialEigenVal, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const int dim = p_space->dim();

  const int numUniqueEigenVals = 3;
  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    V_S(diag.ptr(), ST::one());
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag.ptr())(), null, 0, numUniqueEigenVals, 0 );
    out << "diag =\n" << *diag;
    model->setDiagonalVector(diag);
  }

  const RCP<NonlinearCG<Scalar> > cgSolver =
    createNonlinearCGSolver<Scalar>(model, rcpFromRef(out));

  const RCP<ParameterList> pl = cgSolver->getNonconstParameterList();
  const int minIters = numUniqueEigenVals;
  const int maxIters = minIters + 1;
  //pl->set("Min Num Iterations", minIters);
  pl->set("Max Num Iterations", maxIters);
  cgSolver->setParameterList(pl);

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  ScalarMag g_opt = -1.0;
  const ScalarMag tol = as<Scalar>(g_solve_tol_scale * dim) * ST::eps();
  const ScalarMag alpha_init = 10.0;
  int numIters = -1;
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(tol), optInArg(alpha_init), outArg(numIters) );

  out << "\n";
 
  const ScalarMag err_tol = as<Scalar>(g_error_tol_scale * dim) * ST::eps();
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_COMPARE( numIters, <=, maxIters );
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, FR_partialEigenVal )


//
// Test convergence in full iterations for unique eigen values
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_fullEigenVal, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);
  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const int dim = p_space->dim();
  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag.ptr())(), null );
    out << "diag =\n" << *diag;
    model->setDiagonalVector(diag);
  }

  const RCP<NonlinearCG<Scalar> > cgSolver =
    createNonlinearCGSolver<Scalar>(model, rcpFromRef(out));

  const RCP<ParameterList> pl = cgSolver->getNonconstParameterList();
  const int minIters = dim;
  const int maxIters = minIters+2;
  //pl->set("Min Num Iterations", minIters);
  pl->set("Max Num Iterations", maxIters);
  cgSolver->setParameterList(pl);

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  ScalarMag g_opt = -1.0;
  const ScalarMag tol = as<Scalar>(g_solve_tol_scale * dim) * ST::eps();
  const ScalarMag alpha_init = 10.0;
  int numIters = -1;
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(tol), optInArg(alpha_init), outArg(numIters) );

  out << "\n";
 
  const ScalarMag err_tol = as<Scalar>(g_error_tol_scale * dim) * ST::eps();
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_COMPARE( numIters, <=, maxIters );
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, FR_fullEigenVal )


//
// Test convergence for all unique eigen values but using a scalar product
// that has the effect of clustering the eignvalues seen by the nonlinear CG
// ANA.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_fullEigenValScalarProd, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const int dim = p_space->dim();

  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag.ptr())(), null );
    out << "diag =\n" << *diag;
    model->setDiagonalVector(diag);
  }

  const int numUniqueEigenVals = 3;
  {
    const RCP<VectorBase<Scalar> > diag_bar = createMember(p_space);
    V_S(diag_bar.ptr(), ST::one());
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag_bar.ptr())(), null, 0, numUniqueEigenVals, 0 );
    //applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
    //  null, tuple(diag_bar.ptr())(), null );
    out << "diag_bar =\n" << *diag_bar;
    model->setDiagonalBarVector(diag_bar);
  }

  const RCP<NonlinearCG<Scalar> > cgSolver =
    createNonlinearCGSolver<Scalar>(model, rcpFromRef(out));

  const RCP<ParameterList> pl = cgSolver->getNonconstParameterList();
  const int minIters = numUniqueEigenVals;
  const int maxIters = minIters + 2;
  pl->set("Max Num Iterations", maxIters);
  cgSolver->setParameterList(pl);

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  ScalarMag g_opt = -1.0;
  const ScalarMag tol = as<Scalar>(g_solve_tol_scale * dim) * ST::eps();
  const ScalarMag alpha_init = 10.0;
  int numIters = -1;
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(tol), optInArg(alpha_init), outArg(numIters) );

  out << "\n";
 
  const ScalarMag err_tol = as<Scalar>(g_error_tol_scale * dim) * ST::eps();
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_COMPARE( numIters, <=, maxIters );
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG, FR_fullEigenValScalarProd )


//
// Test general convergence for a general nonlinear objective
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_generalNonlinearProblem, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);

  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag.ptr())(), null );
    out << "diag =\n" << *diag;
    model->setDiagonalVector(diag);
  }

  const ScalarMag nonlinearTermFactor = as<ScalarMag>(g_nonlin_term_factor);
  model->setNonlinearTermFactor(nonlinearTermFactor);

  const RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();
  const RCP<ParameterList> lsPL = parameterList();
  lsPL->set("Min Num Iterations", 1); // Force at least one line search iteration!
  linesearch->setParameterList(lsPL);

  const RCP<NonlinearCG<Scalar> > cgSolver =
    nonlinearCG<Scalar>(model, 0, 0, linesearch);

  cgSolver->setOStream(rcpFromRef(out));

  const RCP<ParameterList> pl = parameterList();
  pl->set("Reinitlaize Linesearch Step Length", true);
  cgSolver->setParameterList(pl);

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  ScalarMag g_opt = -1.0;
  const ScalarMag tol = as<ScalarMag>(g_nonlin_solve_tol);
  const ScalarMag alpha_init = 5.0;
  int numIters = -1;
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      optInArg(tol), optInArg(tol), optInArg(alpha_init), outArg(numIters) );

  out << "\n";
 
  const ScalarMag err_tol = as<ScalarMag>(g_nonlin_error_tol);
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG,
  FR_generalNonlinearProblem )


//
// Test general convergence for a general nonlinear objective passing all
// control options through the PL
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NonlinearCG, FR_generalNonlinearProblem_PL, Scalar )
{

  using Teuchos::optInArg;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const ScalarMag g_offset = as<ScalarMag>(5.0);
  const RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> > model =
    createModel<Scalar>(g_globalDim, g_offset);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);

  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    applyOp<Scalar>( TOpAssignValToGlobalIndex<Scalar>(),
      null, tuple(diag.ptr())(), null );
    out << "diag =\n" << *diag;
    model->setDiagonalVector(diag);
  }

  const ScalarMag nonlinearTermFactor = as<ScalarMag>(g_nonlin_term_factor);
  model->setNonlinearTermFactor(nonlinearTermFactor);

  const RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();
  const RCP<ParameterList> lsPL = parameterList();
  lsPL->set("Min Num Iterations", 1); // Force at least one line search iteration!
  linesearch->setParameterList(lsPL);

  const RCP<NonlinearCG<Scalar> > cgSolver =
    nonlinearCG<Scalar>(model, 0, 0, linesearch);

  cgSolver->setOStream(rcpFromRef(out));

  const RCP<VectorBase<Scalar> > p = createMember(p_space);
  V_S( p.ptr(), ST::zero() );
  
  const double tol = as<double>(g_nonlin_solve_tol);
  const double alpha_init = as<double>(5.0);
  const RCP<ParameterList> pl = parameterList();
  pl->set("Initial Linesearch Step Length", alpha_init);
  pl->set("Reinitlaize Linesearch Step Length", true);
  pl->set("Objective Reduction Tol", tol);
  pl->set("Objective Gradient Tol", tol);
  cgSolver->setParameterList(pl);

  ScalarMag g_opt = -1.0;
  int numIters = -1;
  const NCGU::ESolveReturn solveResult =
    cgSolver->doSolve( p.ptr(), outArg(g_opt),
      null, null, null, outArg(numIters) );
  
  out << "\n";
 
  const ScalarMag err_tol = as<ScalarMag>(g_nonlin_error_tol);
  TEST_EQUALITY(solveResult, NCGU::SOLVE_SOLUTION_FOUND);
  TEST_FLOATING_EQUALITY(g_opt, g_offset, err_tol);
  const bool result = Thyra::testRelNormDiffErr<Scalar>(
    "p", *p,
    "ps", *model->getSolutionVector(),
    "err_tol", err_tol,
    "2*err_tol", as<ScalarMag>(2.0)*err_tol,
    &out
    );
  if (!result) success = false;
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( NonlinearCG,
  FR_generalNonlinearProblem_PL )


} // namespace


