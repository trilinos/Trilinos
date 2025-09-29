// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DiagonalQuadraticResponseOnlyModelEvaluator.hpp"
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

int g_localDim = 4;

double g_tol_scale = 10.0;

TEUCHOS_STATIC_SETUP() {
  Teuchos::UnitTestRepository::getCLP().setOption(
      "local-dim", &g_localDim, "Number of local vector elements on each process");
  Teuchos::UnitTestRepository::getCLP().setOption(
      "tol-scale", &g_tol_scale, "Floating point tolerance scaling");
}

template <class Scalar>
inline Scalar sqr(const Scalar &x) { return x * x; }

template <class Scalar>
inline Scalar cube(const Scalar &x) { return x * x * x; }

//
// Unit tests
//

//
// Test basic construction
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DiagonalQuadraticResponseOnlyModelEvaluator,
                                  basic, Scalar) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  // typedef Thyra::ModelEvaluatorBase MEB; // unused
  using Thyra::create_DgDp_mv;
  using Thyra::derivativeGradient;
  using Thyra::eval_g_DgDp;
  using Thyra::get_ele;
  using Thyra::get_mv;
  using Thyra::norm_2;

  RCP<const Thyra::ModelEvaluator<Scalar> >
      model = Thyra::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);

  TEST_ASSERT(!is_null(model));
  TEST_EQUALITY_CONST(model->Np(), 1);
  TEST_EQUALITY_CONST(model->Ng(), 1);

  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  RCP<const Thyra::VectorSpaceBase<Scalar> > g_space = model->get_g_space(0);

  RCP<Thyra::VectorBase<Scalar> > p_init = createMember(p_space);
  const Scalar val                       = as<Scalar>(2.0);
  out << "\nval = " << val << "\n";
  Thyra::V_S(p_init.ptr(), val);

  RCP<Thyra::VectorBase<Scalar> >
      g      = createMember(g_space),
      g_grad = createMember(p_space);

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
                      g.ptr(), derivativeGradient<Scalar>(g_grad));

  out << "\ng =\n"
      << *g;
  out << "\ng_grad =\n"
      << *g_grad;

  const Ordinal globalDim = p_space->dim();
  out << "\nglobalDim = " << globalDim << "\n";

  TEST_FLOATING_EQUALITY(get_ele<Scalar>(*g, 0),
                         as<Scalar>(0.5 * globalDim) * val * val, as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));

  TEST_FLOATING_EQUALITY(
      norm_2<Scalar>(*g_grad),
      ST::magnitude(ST::squareroot(as<Scalar>(globalDim) * val * val)),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
    DiagonalQuadraticResponseOnlyModelEvaluator, basic)

//
// Test with all of knobs setup
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DiagonalQuadraticResponseOnlyModelEvaluator,
                                  offsets, Scalar) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  // typedef Thyra::ModelEvaluatorBase MEB; // unused
  using Thyra::derivativeGradient;
  using Thyra::eval_g_DgDp;
  using Thyra::get_ele;
  using Thyra::get_mv;
  using Thyra::norm_2;
  using Thyra::V_S;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;

  const RCP<Thyra::DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
      model = Thyra::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<const VectorSpaceBase<Scalar> > g_space = model->get_g_space(0);

  const Scalar p_soln_val = as<Scalar>(3.0);
  {
    const RCP<VectorBase<Scalar> > p_soln = createMember(p_space);
    V_S(p_soln.ptr(), p_soln_val);
    model->setSolutionVector(p_soln);
  }

  const Scalar diag_val = as<Scalar>(4.0);
  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    V_S(diag.ptr(), diag_val);
    model->setDiagonalVector(diag);
  }

  const Scalar g_offset = as<Scalar>(5.0);
  {
    model->setScalarOffset(g_offset);
  }

  const Scalar nonlinearTermFactor = 1e-3;
  {
    model->setNonlinearTermFactor(nonlinearTermFactor);
  }

  const Scalar p_val                    = as<Scalar>(2.0);
  const RCP<VectorBase<Scalar> > p_init = createMember(p_space);
  V_S(p_init.ptr(), p_val);

  RCP<VectorBase<Scalar> >
      g      = createMember(g_space),
      g_grad = createMember(p_space);

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
                      g.ptr(), derivativeGradient<Scalar>(g_grad));

  out << "\ng =\n"
      << *g;
  out << "\ng_grad =\n"
      << *g_grad;

  const Ordinal globalDim = p_space->dim();
  out << "\nglobalDim = " << globalDim << "\n";

  TEST_FLOATING_EQUALITY(
      get_ele<Scalar>(*g, 0),
      as<Scalar>(
          0.5 * globalDim * (diag_val * sqr(p_val - p_soln_val) + nonlinearTermFactor * cube(p_val - p_soln_val)) + g_offset),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));

  TEST_FLOATING_EQUALITY(
      sum(*g_grad),
      as<Scalar>(
          globalDim * (diag_val * (p_val - p_soln_val) + 1.5 * nonlinearTermFactor * sqr(p_val - p_soln_val))),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
    DiagonalQuadraticResponseOnlyModelEvaluator, offsets)

//
// Test non-Euclidean scalar product
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DiagonalQuadraticResponseOnlyModelEvaluator,
                                  nonEuclideanScalarProd, Scalar) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Ordinal;
  // typedef Thyra::ModelEvaluatorBase MEB; // unused
  using Thyra::derivativeGradient;
  using Thyra::eval_g_DgDp;
  using Thyra::get_ele;
  using Thyra::get_mv;
  using Thyra::norm_2;
  using Thyra::scalarProd;
  using Thyra::V_S;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;

  const RCP<Thyra::DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
      model = Thyra::diagonalQuadraticResponseOnlyModelEvaluator<Scalar>(g_localDim);

  const RCP<const VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<const VectorSpaceBase<Scalar> > g_space = model->get_g_space(0);

  const Scalar diag_val = as<Scalar>(4.0);
  {
    const RCP<VectorBase<Scalar> > diag = createMember(p_space);
    V_S(diag.ptr(), diag_val);
    model->setDiagonalVector(diag);
  }

  const Scalar diag_bar_val = as<Scalar>(3.0);
  {
    const RCP<VectorBase<Scalar> > diag_bar = createMember(p_space);
    V_S(diag_bar.ptr(), diag_bar_val);
    model->setDiagonalBarVector(diag_bar);
  }

  const Scalar g_offset = as<Scalar>(5.0);
  {
    model->setScalarOffset(g_offset);
  }

  const Scalar nonlinearTermFactor = 1e-3;
  {
    model->setNonlinearTermFactor(nonlinearTermFactor);
  }

  const Scalar x_val               = as<Scalar>(5.0);
  const RCP<VectorBase<Scalar> > x = createMember(p_space);
  V_S(x.ptr(), x_val);

  const Scalar y_val               = as<Scalar>(6.0);
  const RCP<VectorBase<Scalar> > y = createMember(p_space);
  V_S(y.ptr(), y_val);

  const Ordinal globalDim = p_space->dim();
  out << "\nglobalDim = " << globalDim << "\n";

  TEST_FLOATING_EQUALITY(
      scalarProd(*x, *y),
      as<Scalar>(globalDim * (diag_val / diag_bar_val) * x_val * y_val),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));

  const Scalar p_soln_val = as<Scalar>(3.0);
  {
    const RCP<VectorBase<Scalar> > p_soln = createMember(p_space);
    V_S(p_soln.ptr(), p_soln_val);
    model->setSolutionVector(p_soln);
  }

  const Scalar p_val                    = as<Scalar>(2.0);
  const RCP<VectorBase<Scalar> > p_init = createMember(p_space);
  V_S(p_init.ptr(), p_val);

  RCP<VectorBase<Scalar> >
      g      = createMember(g_space),
      g_grad = createMember(p_space);

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
                      g.ptr(), derivativeGradient<Scalar>(g_grad));

  out << "\ng =\n"
      << *g;
  out << "\ng_grad =\n"
      << *g_grad;

  eval_g_DgDp<Scalar>(*model, 0, *p_init, 0,
                      g.ptr(), derivativeGradient<Scalar>(g_grad));

  TEST_FLOATING_EQUALITY(
      get_ele<Scalar>(*g, 0),
      as<Scalar>(
          0.5 * globalDim * (diag_val * sqr(p_val - p_soln_val) + nonlinearTermFactor * cube(p_val - p_soln_val)) + g_offset),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));

  TEST_FLOATING_EQUALITY(
      sum(*g_grad),
      as<Scalar>(
          globalDim * (diag_bar_val / diag_val) * (diag_val * (p_val - p_soln_val) + 1.5 * nonlinearTermFactor * sqr(p_val - p_soln_val))),
      as<ScalarMag>(g_tol_scale * ST::eps() / globalDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
    DiagonalQuadraticResponseOnlyModelEvaluator, nonEuclideanScalarProd)

}  // namespace
