//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_DAEStructure_hpp
#define Tempus_DAEStructure_hpp

#include <algorithm>
#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

/** \brief Structure of an index-1 differential-algebraic system.
 *
 *  The equations f(x_dot, x, p, t) = 0 are split into differential and
 *  algebraic equations.  The algebraic equations are rows of the mass matrix
 *  E = df/dx_dot that are zero, e.g. constraints solved at each time
 *  (the velocity equations of a coupled velocity-thickness ice sheet model),
 *  or Dirichlet conditions written as f_i = x_i - d_i(p,t).  The structure is
 *  given by a mask in the residual space, equal to 1 on the algebraic
 *  equations and 0 elsewhere: only the equations are flagged, not the
 *  unknowns.  With A = df/dx and the diagonal masks P_a and P_d = I - P_a of
 *  the algebraic and differential equations, the system has index 1 when the
 *  matrix
 *      C = P_a A + P_d E
 *  (rows of df/dx for the algebraic equations, rows of df/dx_dot for the
 *  others) is nonsingular: differentiating the algebraic equations once in
 *  time gives C x_dot = -(P_a df/dt + P_d f(0, x, p, t)), which determines
 *  x_dot.  Semi-explicit DAEs, x = (x_d, x_a) with
 *  E = [E_dd 0; 0 0] and df_a/dx_a nonsingular, are the typical case.
 *
 *  The structure is given to Tempus with the parameter "Algebraic Equations
 *  Mask" of the "Sensitivities" sublist, a
 *  Teuchos::RCP<const Thyra::VectorBase<Scalar> > (see getDAEStructure()).
 *  It is currently used by IntegratorAdjointSensitivity, to compute the
 *  jumps of the adjoint of "Final Time" and "Sum Over Times" responses (see
 *  AdjointDAEJumpSolver).  The forward integration does not use it: the
 *  forward stepper must be suited to index-1 DAEs (e.g. Backward Euler) and
 *  must not compute consistent initial conditions.
 */
template <class Scalar>
class DAEStructure {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;
  typedef Thyra::MultiVectorBase<Scalar> MultiVector;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

  /// Name of the parameter giving the mask of the algebraic equations
  static const std::string& algebraicMaskName()
  {
    static const std::string name("Algebraic Equations Mask");
    return name;
  }

  /// Constructor from the mask of the algebraic equations (entries 0 or 1)
  explicit DAEStructure(const Teuchos::RCP<const Vector>& algebraic_mask)
    : algebraic_mask_(algebraic_mask)
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEUCHOS_TEST_FOR_EXCEPTION(
        algebraic_mask_ == Teuchos::null, std::invalid_argument,
        "Error - Tempus::DAEStructure: null \"" << algebraicMaskName() << "\".\n");

    // The entries must be 0 or 1: m.*m - m = 0
    Teuchos::RCP<Vector> w = algebraic_mask_->clone_v();
    Thyra::ele_wise_scale(*algebraic_mask_, w.ptr());
    Thyra::Vp_StV(w.ptr(), -ST::one(), *algebraic_mask_);
    TEUCHOS_TEST_FOR_EXCEPTION(
        Thyra::norm_inf(*w) != Teuchos::ScalarTraits<Mag>::zero(),
        std::invalid_argument,
        "Error - Tempus::DAEStructure: the entries of \""
            << algebraicMaskName()
            << "\" must be 1 (algebraic equations) or 0 (other equations).\n");

    Teuchos::RCP<Vector> d = Thyra::createMember(algebraic_mask_->space());
    Thyra::put_scalar(ST::one(), d.ptr());
    Thyra::Vp_StV(d.ptr(), -ST::one(), *algebraic_mask_);
    differential_mask_ = d;
  }

  /// Mask of the algebraic equations (1 on the algebraic equations)
  Teuchos::RCP<const Vector> getAlgebraicMask() const
  {
    return algebraic_mask_;
  }

  /// Mask of the other equations (1 - mask of the algebraic equations)
  Teuchos::RCP<const Vector> getDifferentialMask() const
  {
    return differential_mask_;
  }

  /// v <- P_d v: zero the algebraic components of each column of v
  void zeroAlgebraicComponents(MultiVector& v) const
  {
    for (int j = 0; j < v.domain()->dim(); ++j)
      Thyra::ele_wise_scale(*differential_mask_, v.col(j).ptr());
  }

  /// v <- P_a v: zero the differential components of each column of v
  void zeroDifferentialComponents(MultiVector& v) const
  {
    for (int j = 0; j < v.domain()->dim(); ++j)
      Thyra::ele_wise_scale(*algebraic_mask_, v.col(j).ptr());
  }

  /** \brief Check that the rows of the algebraic equations of the mass matrix
   *  df/dx_dot (model evaluated with alpha = 1, beta = 0 at inArgs) are zero.
   *
   *  The mass matrix is applied to a random vector.  Models without W_op are
   *  not checked.
   */
  void checkMassMatrix(const Thyra::ModelEvaluator<Scalar>& model,
                       const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                       const std::string& caller) const
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::OutArgs<Scalar> outArgs = model.createOutArgs();
    if (!outArgs.supports(MEB::OUT_ARG_W_op)) return;

    MEB::InArgs<Scalar> massInArgs = inArgs;
    if (massInArgs.supports(MEB::IN_ARG_alpha))
      massInArgs.set_alpha(Scalar(1.0));
    if (massInArgs.supports(MEB::IN_ARG_beta)) massInArgs.set_beta(Scalar(0.0));
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > W_op = model.create_W_op();
    outArgs.set_W_op(W_op);
    model.evalModel(massInArgs, outArgs);

    Teuchos::RCP<Vector> v = Thyra::createMember(model.get_x_space());
    Thyra::randomize(Scalar(0.5), Scalar(1.5), v.ptr());
    Teuchos::RCP<Vector> r = Thyra::createMember(model.get_f_space());
    W_op->apply(Thyra::NOTRANS, *v, r.ptr(), Scalar(1.0), Scalar(0.0));
    const Mag scale = Thyra::norm_inf(*r);
    Thyra::ele_wise_scale(*algebraic_mask_, r.ptr());  // algebraic rows only
    const Mag err = Thyra::norm_inf(*r);
    TEUCHOS_TEST_FOR_EXCEPTION(
        err > Mag(1.0e-12) * scale, std::logic_error,
        "Error - " << caller << ": the rows of the mass matrix df/dx_dot "
                   << "(model evaluated with alpha = 1, beta = 0) flagged by \""
                   << algebraicMaskName() << "\" are not zero (|E v| = " << err
                   << " on these rows, " << scale
                   << " on all rows).  Only equations that do not depend on "
                   << "x_dot can be flagged as algebraic.\n");
  }

 private:
  Teuchos::RCP<const Vector> algebraic_mask_;
  Teuchos::RCP<const Vector> differential_mask_;
};

/** \brief The DAE structure given by the "Algebraic Equations Mask"
 *  parameter of pl, or null if the parameter is not set (or null).
 *
 *  Both Teuchos::RCP<const Thyra::VectorBase<Scalar> > and
 *  Teuchos::RCP<Thyra::VectorBase<Scalar> > are accepted.
 */
template <class Scalar>
Teuchos::RCP<const DAEStructure<Scalar> > getDAEStructure(
    const Teuchos::ParameterList& pl)
{
  typedef Teuchos::RCP<const Thyra::VectorBase<Scalar> > ConstVec;
  typedef Teuchos::RCP<Thyra::VectorBase<Scalar> > Vec;
  const std::string& name = DAEStructure<Scalar>::algebraicMaskName();
  if (!pl.isParameter(name)) return Teuchos::null;
  ConstVec mask;
  if (pl.isType<ConstVec>(name))
    mask = pl.get<ConstVec>(name);
  else if (pl.isType<Vec>(name))
    mask = pl.get<Vec>(name);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        "Error - the parameter \""
            << name
            << "\" must be a Teuchos::RCP<const Thyra::VectorBase<Scalar> >.\n");
  if (mask == Teuchos::null) return Teuchos::null;
  return Teuchos::rcp(new DAEStructure<Scalar>(mask));
}

}  // namespace Tempus

#endif  // Tempus_DAEStructure_hpp
