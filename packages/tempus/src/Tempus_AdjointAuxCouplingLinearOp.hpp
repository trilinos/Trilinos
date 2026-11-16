//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AdjointAuxCouplingLinearOp_hpp
#define Tempus_AdjointAuxCouplingLinearOp_hpp

#include "Teuchos_Assert.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

/** \brief Coupling block s*op(df/dp) of the W of
 *  AdjointAuxSensitivityModelEvaluator.
 *
 *  The auxiliary equations dz/dt - df/dp^T*y = 0 of the adjoint sensitivity
 *  equations depend on the adjoint y through -df/dp^T*y, so the block of W of
 *  the auxiliary rows and the adjoint columns is -beta*df/dp^T.  This
 *  operator is s*op(df/dp), with df/dp evaluated by the model in one of the
 *  forms of applyDfDpTranspose() (see AdjointSensitivityUtils::evalDfDp()):
 *  a linear operator or a Jacobian multivector, applied with
 *  op = CONJTRANS, or a gradient multivector, applied with op = CONJ.  It
 *  maps the residual (f) space of the model to the space of the parameter.
 *  Before initialize(), and with a null df/dp, the operator is zero.
 *
 *  The left and right scalings by diagonal matrices are supported, so that
 *  the operator can be a block of an operator scaled by a nonlinear solver
 *  (e.g., NOX with "Function Scaling" = "Row Sum"); the scalings are stored
 *  with the operator, and initialize() resets them.
 *
 *  The row statistics (absolute row and column sums) are reported as zero:
 *  the block is left out of the row statistics of W, which are those of its
 *  block diagonal.  The row scaling of the auxiliary rows (e.g., the inverse
 *  row sums) is then uniform, as needed by the scaled identity of the
 *  diagonal block of these rows (ScaledIdentityLinearOpWithSolve).
 */
template <class Scalar>
class AdjointAuxCouplingLinearOp
  : virtual public Thyra::RowStatLinearOpBase<Scalar>,
    virtual public Thyra::ScaledLinearOpBase<Scalar> {
 public:
  typedef Teuchos::ScalarTraits<Scalar> ST;

  /** \brief Zero operator from domain (the residual space of the model) to
   *  range (the parameter space). */
  AdjointAuxCouplingLinearOp(
      const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >& range,
      const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >& domain)
    : range_(range), domain_(domain), dfdp_trans_(Thyra::CONJTRANS), s_(ST::zero())
  {
    TEUCHOS_ASSERT(Teuchos::nonnull(range_));
    TEUCHOS_ASSERT(Teuchos::nonnull(domain_));
  }

  /** \brief Set the operator to s*op(dfdp), with op = dfdp_trans, and remove
   *  the scalings.
   *
   *  The operator keeps a pointer to dfdp, it does not copy it.
   */
  void initialize(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >& dfdp,
                  const Thyra::EOpTransp dfdp_trans, const Scalar& s)
  {
    dfdp_          = dfdp;
    dfdp_trans_    = dfdp_trans;
    s_             = s;
    left_scaling_  = Teuchos::null;
    right_scaling_ = Teuchos::null;
  }

  /** @name Overridden from LinearOpBase */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > range() const
  {
    return range_;
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const
  {
    return domain_;
  }

  //@}

 protected:
  /** @name Overridden from LinearOpBase */
  //@{

  bool opSupportedImpl(Thyra::EOpTransp M_trans) const
  {
    if (Teuchos::is_null(dfdp_)) return true;
    // A complex conjugation of the scalings is not implemented
    if (ST::isComplex && isConjugate(M_trans) && isScaled()) return false;
    return Thyra::opSupported(*dfdp_, Thyra::trans_trans(dfdp_trans_, M_trans));
  }

  /** \brief Y = alpha*op_M(D_l*s*op(df/dp)*D_r)*X + beta*Y, with D_l and D_r
   *  the left and right scalings (identity if not set). */
  void applyImpl(const Thyra::EOpTransp M_trans,
                 const Thyra::MultiVectorBase<Scalar>& X,
                 const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> >& Y,
                 const Scalar alpha, const Scalar beta) const
  {
    using Teuchos::RCP;
    typedef Thyra::MultiVectorBase<Scalar> MV;

    if (Teuchos::is_null(dfdp_)) {
      scaleOrZero(beta, *Y);
      return;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        ST::isComplex && isConjugate(M_trans) && isScaled(), std::logic_error,
        "Error - Tempus::AdjointAuxCouplingLinearOp::apply(): a conjugate "
        "apply of the scaled operator is not supported.\n");

    // op_M(D_l*s*op(A)*D_r) = D_l*s*op_M(op(A))*D_r for NOTRANS and CONJ,
    // and D_r*s*op_M(op(A))*D_l for TRANS and CONJTRANS (real scalings)
    const bool notrans = (Thyra::real_trans(M_trans) == Thyra::NOTRANS);
    const RCP<const Thyra::VectorBase<Scalar> > pre_scaling =
        notrans ? right_scaling_ : left_scaling_;
    const RCP<const Thyra::VectorBase<Scalar> > post_scaling =
        notrans ? left_scaling_ : right_scaling_;
    const Scalar s               = isConjugate(M_trans) ? ST::conjugate(s_) : s_;
    const Thyra::EOpTransp trans = Thyra::trans_trans(dfdp_trans_, M_trans);

    RCP<const MV> Xs = Teuchos::rcpFromRef(X);
    if (Teuchos::nonnull(pre_scaling)) {
      RCP<MV> X_scaled = X.clone_mv();
      X_scaled->scaleLeft(*pre_scaling);
      Xs = X_scaled;
    }
    if (Teuchos::is_null(post_scaling)) {
      Thyra::apply(*dfdp_, trans, *Xs, Y, alpha * s, beta);
    }
    else {
      RCP<MV> T = Thyra::createMembers(notrans ? range_ : domain_,
                                       X.domain()->dim());
      Thyra::apply(*dfdp_, trans, *Xs, T.ptr(), s, ST::zero());
      T->scaleLeft(*post_scaling);
      scaleOrZero(beta, *Y);
      Thyra::update(alpha, *T, Y);
    }
  }

  //@}

  /** @name Overridden from RowStatLinearOpBase */
  //@{

  /** \brief The absolute row and column sums are reported as zero (see the
   *  class description), so their inverses are not supported. */
  bool rowStatIsSupportedImpl(
      const Thyra::RowStatLinearOpBaseUtils::ERowStat rowStat) const
  {
    return rowStat == Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM ||
           rowStat == Thyra::RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
  }

  void getRowStatImpl(
      const Thyra::RowStatLinearOpBaseUtils::ERowStat rowStat,
      const Teuchos::Ptr<Thyra::VectorBase<Scalar> >& rowStatVec) const
  {
    TEUCHOS_ASSERT(this->rowStatIsSupported(rowStat));
    Thyra::put_scalar(ST::zero(), rowStatVec);
  }

  //@}

  /** @name Overridden from ScaledLinearOpBase */
  //@{

  bool supportsScaleLeftImpl() const { return true; }

  bool supportsScaleRightImpl() const { return true; }

  void scaleLeftImpl(const Thyra::VectorBase<Scalar>& row_scaling)
  {
    accumulateScaling(row_scaling, left_scaling_);
  }

  void scaleRightImpl(const Thyra::VectorBase<Scalar>& col_scaling)
  {
    accumulateScaling(col_scaling, right_scaling_);
  }

  //@}

 private:
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > range_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > domain_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > dfdp_;
  Thyra::EOpTransp dfdp_trans_;
  Scalar s_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > left_scaling_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > right_scaling_;

  bool isScaled() const
  {
    return Teuchos::nonnull(left_scaling_) || Teuchos::nonnull(right_scaling_);
  }

  static bool isConjugate(const Thyra::EOpTransp M_trans)
  {
    return M_trans == Thyra::CONJ || M_trans == Thyra::CONJTRANS;
  }

  /// Y = beta*Y, with Y = 0 (whatever its entries) if beta = 0
  static void scaleOrZero(const Scalar beta, Thyra::MultiVectorBase<Scalar>& Y)
  {
    if (beta == ST::zero())
      Thyra::assign(Teuchos::ptrFromRef(Y), ST::zero());
    else if (beta != ST::one())
      Thyra::scale(beta, Teuchos::ptrFromRef(Y));
  }

  /// scaling <- diag(d)*scaling (scaling = d if not set)
  static void accumulateScaling(const Thyra::VectorBase<Scalar>& d,
                                Teuchos::RCP<Thyra::VectorBase<Scalar> >& scaling)
  {
    if (Teuchos::is_null(scaling))
      scaling = d.clone_v();
    else
      Thyra::ele_wise_scale(d, scaling.ptr());
  }
};

}  // namespace Tempus

#endif  // Tempus_AdjointAuxCouplingLinearOp_hpp
