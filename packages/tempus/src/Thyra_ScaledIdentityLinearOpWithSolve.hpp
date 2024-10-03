//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_ScaledIdentityLinearOpWithSolve_hpp
#define Thyra_ScaledIdentityLinearOpWithSolve_hpp

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Thyra {

/**
 * \brief Implicit concrete <tt>LinearOpBase</tt> subclass that
 * takes a flattended out multi-vector and performs a multi-RHS apply with it.
 */
template <class Scalar>
class ScaledIdentityLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  ScaledIdentityLinearOpWithSolve() {}

  void initialize(const RCP<const VectorSpaceBase<Scalar> >& space,
                  const Scalar& s)
  {
    validateInitialize(space);
    space_ = space;
    s_     = s;
  }

  void uninitialize() { space_ = Teuchos::null; }

  RCP<const VectorSpaceBase<Scalar> > space() const { return space_; }
  Scalar scale() const { return s_; }
  void setScale(const Scalar& s) { s_ = s; }

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  RCP<const VectorSpaceBase<Scalar> > range() const { return space_; }

  RCP<const VectorSpaceBase<Scalar> > domain() const { return space_; }

  RCP<const LinearOpBase<Scalar> > clone() const
  {
    RCP<ScaledIdentityLinearOpWithSolve<Scalar> > op =
        rcp(new ScaledIdentityLinearOpWithSolve<Scalar>());
    op->initialize(space_, s_);
    return op;
  }
  //@}

 protected:
  /** @name Overridden from LinearOpBase */
  //@{
  bool opSupportedImpl(EOpTransp /* M_trans */) const { return true; }

  void applyImpl(const EOpTransp M_trans, const MultiVectorBase<Scalar>& X,
                 const Ptr<MultiVectorBase<Scalar> >& Y, const Scalar alpha,
                 const Scalar beta) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::scale(beta, Y);
    if (M_trans == CONJ || M_trans == CONJTRANS)
      V_StVpV(Y, ST::conjugate(s_) * alpha, X, *Y);
    else
      V_StVpV(Y, s_ * alpha, X, *Y);
  }
  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  bool solveSupportsImpl(EOpTransp /* M_trans */) const { return true; }

  bool solveSupportsNewImpl(
      EOpTransp /* M_trans */,
      const Ptr<const SolveCriteria<Scalar> > /* solveCriteria */) const
  {
    return true;
  }

  bool solveSupportsSolveMeasureTypeImpl(
      EOpTransp /* M_trans */,
      const SolveMeasureType& /* solveMeasureType */) const
  {
    return true;
  }

  SolveStatus<Scalar> solveImpl(
      const EOpTransp M_trans, const MultiVectorBase<Scalar>& B,
      const Ptr<MultiVectorBase<Scalar> >& X,
      const Ptr<const SolveCriteria<Scalar> > /* solveCriteria */) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    assign(X, ST::zero());
    if (M_trans == CONJ || M_trans == CONJTRANS)
      V_StVpV(X, ST::one() / ST::conjugate(s_), B, *X);
    else
      V_StVpV(X, ST::one() / s_, B, *X);
    SolveStatus<Scalar> solveStatus;
    solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
    return solveStatus;
  }
  //@}

 private:
  // //////////////////////////////
  // Private data members

  RCP<const VectorSpaceBase<Scalar> > space_;
  Scalar s_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const VectorSpaceBase<Scalar> >& space)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(space));
#else
    (void)space;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates ScaledIdentityLinearOpWithSolve
 */
template <class Scalar>
RCP<ScaledIdentityLinearOpWithSolve<Scalar> > scaledIdentity()
{
  return Teuchos::rcp(new ScaledIdentityLinearOpWithSolve<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates ScaledIdentityLinearOpWithSolve
 */
template <class Scalar>
RCP<ScaledIdentityLinearOpWithSolve<Scalar> > scaledIdentity(
    const RCP<const VectorSpaceBase<Scalar> >& space, const Scalar& s)
{
  RCP<ScaledIdentityLinearOpWithSolve<Scalar> > op =
      Teuchos::rcp(new ScaledIdentityLinearOpWithSolve<Scalar>());
  op->initialize(space, s);
  return op;
}

}  // end namespace Thyra

#endif
