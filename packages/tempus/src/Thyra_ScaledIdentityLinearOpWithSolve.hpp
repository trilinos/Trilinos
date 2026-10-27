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

#include <string>

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Thyra {

/**
 * \brief Scaled identity operator s*I, with its solve.
 *
 * The row statistics (absolute row and column sums |s|) and the left and
 * right scaling by uniform vectors (all entries equal to d, s <- d*s) are
 * supported, so that the operator can be a block of an operator scaled by a
 * nonlinear solver, e.g., NOX with "Function Scaling" = "Row Sum" (the
 * auxiliary block of AdjointAuxSensitivityModelEvaluator).  A scaling by a
 * nonuniform vector throws, since the result is not a scaled identity.
 */
template <class Scalar>
class ScaledIdentityLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar>,
    virtual public RowStatLinearOpBase<Scalar>,
    virtual public ScaledLinearOpBase<Scalar> {
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

  /** @name Overridden from RowStatLinearOpBase */
  //@{

  /** \brief The absolute row and column sums are |s|; their inverses
   *  require s != 0. */
  bool rowStatIsSupportedImpl(
      const RowStatLinearOpBaseUtils::ERowStat rowStat) const
  {
    switch (rowStat) {
      case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM: return true;
      case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
      case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
        return s_ != Teuchos::ScalarTraits<Scalar>::zero();
    }
    return false;
  }

  void getRowStatImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat,
                      const Ptr<VectorBase<Scalar> >& rowStatVec) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEUCHOS_ASSERT(this->rowStatIsSupported(rowStat));
    const Scalar abs_s = ST::magnitude(s_);
    if (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM ||
        rowStat == RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM)
      put_scalar(ST::one() / abs_s, rowStatVec);
    else
      put_scalar(abs_s, rowStatVec);
  }
  //@}

  /** @name Overridden from ScaledLinearOpBase */
  //@{

  bool supportsScaleLeftImpl() const { return true; }

  bool supportsScaleRightImpl() const { return true; }

  /** \brief diag(d) * s*I = (d*s) I, for a uniform vector d. */
  void scaleLeftImpl(const VectorBase<Scalar>& row_scaling)
  {
    s_ *= uniformValue(row_scaling, "scaleLeft");
  }

  /** \brief s*I * diag(d) = (d*s) I, for a uniform vector d. */
  void scaleRightImpl(const VectorBase<Scalar>& col_scaling)
  {
    s_ *= uniformValue(col_scaling, "scaleRight");
  }
  //@}

 private:
  // //////////////////////////////
  // Private data members

  RCP<const VectorSpaceBase<Scalar> > space_;
  Scalar s_;

  // //////////////////////////////
  // Private member functions

  /// The value d of the entries of a uniform vector (throws otherwise)
  static Scalar uniformValue(const VectorBase<Scalar>& v,
                             const std::string& caller)
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    if (v.space()->dim() == 0) return ST::one();
    const Scalar d                = get_ele(v, 0);
    RCP<VectorBase<Scalar> > diff = v.clone_v();
    add_scalar(-d, diff.ptr());
    TEUCHOS_TEST_FOR_EXCEPTION(
        norm_inf(*diff) > Mag(100) * ST::eps() * ST::magnitude(d),
        std::logic_error,
        "Error - Thyra::ScaledIdentityLinearOpWithSolve::"
            << caller
            << ": only a scaling by a uniform vector (all entries equal) "
               "keeps the operator a scaled identity.\n");
    return d;
  }

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
