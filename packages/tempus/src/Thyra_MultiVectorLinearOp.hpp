//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_MultiVectorLinearOp_hpp
#define Thyra_MultiVectorLinearOp_hpp

#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

/**
 * \brief Implicit concrete <tt>LinearOpBase</tt> subclass that
 * takes a flattended out multi-vector and performs a multi-RHS apply with it.
 */
template <class Scalar>
class MultiVectorLinearOp : virtual public RowStatLinearOpBase<Scalar>,
                            virtual public ScaledLinearOpBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  MultiVectorLinearOp() {}

  void nonconstInitialize(
      const RCP<LinearOpBase<Scalar> > &op,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
    validateInitialize(op, multiVecRange, multiVecDomain);
    op_             = op;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  void initialize(const RCP<const LinearOpBase<Scalar> > &op,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecRange,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecDomain)
  {
    validateInitialize(op, multiVecRange, multiVecDomain);
    op_             = op;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  RCP<LinearOpBase<Scalar> > getNonconstLinearOp()
  {
    return op_.getNonconstObj();
  }

  RCP<const LinearOpBase<Scalar> > getLinearOp() const
  {
    return op_.getConstObj();
  }

  void uninitialize()
  {
    op_.uninitialize();
    multiVecRange_  = Teuchos::null;
    multiVecDomain_ = Teuchos::null;
  }

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  RCP<const VectorSpaceBase<Scalar> > range() const { return multiVecRange_; }

  RCP<const VectorSpaceBase<Scalar> > domain() const { return multiVecDomain_; }

  RCP<const LinearOpBase<Scalar> > clone() const
  {
    return Teuchos::null;  // ToDo: Implement if needed ???
  }
  //@}

 protected:
  /** @name Overridden from LinearOpBase */
  //@{
  bool opSupportedImpl(EOpTransp M_trans) const
  {
    return Thyra::opSupported(*op_.getConstObj(), M_trans);
  }

  void applyImpl(const EOpTransp M_trans, const MultiVectorBase<Scalar> &XX,
                 const Ptr<MultiVectorBase<Scalar> > &YY, const Scalar alpha,
                 const Scalar beta) const
  {
    using Teuchos::dyn_cast;
    typedef DefaultMultiVectorProductVector<Scalar> MVPV;

    const Ordinal numCols = XX.domain()->dim();

    for (Ordinal col_j = 0; col_j < numCols; ++col_j) {
      const RCP<const VectorBase<Scalar> > x = XX.col(col_j);
      const RCP<VectorBase<Scalar> > y       = YY->col(col_j);

      RCP<const MultiVectorBase<Scalar> > X =
          dyn_cast<const MVPV>(*x).getMultiVector().assert_not_null();
      RCP<MultiVectorBase<Scalar> > Y =
          dyn_cast<MVPV>(*y).getNonconstMultiVector().assert_not_null();

      Thyra::apply(*op_.getConstObj(), M_trans, *X, Y.ptr(), alpha, beta);
    }
  }
  //@}

  /** @name Overridden from RowStatLinearOpBase */
  //@{

  /** \brief Determine if a given row stat is supported. */
  bool rowStatIsSupportedImpl(
      const RowStatLinearOpBaseUtils::ERowStat rowStat) const
  {
    using Teuchos::rcp_dynamic_cast;
    return rcp_dynamic_cast<const RowStatLinearOpBase<Scalar> >(
               op_.getConstObj())
        ->rowStatIsSupported(rowStat);
  }

  /** \brief Get some statistics about a supported row.
   *
   * \pre <tt>this->rowStatIsSupported(rowStat)==true</tt>
   */
  void getRowStatImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat,
                      const Ptr<VectorBase<Scalar> > &rowStatVec) const
  {
    TEUCHOS_ASSERT(this->rowStatIsSupported(rowStat));

    // Compute the scaling vector from the underlying operator and assign
    // it to each column of the scaling multivector.  We only use the first
    // column in the scaling below, but this makes the scaling vector
    // consistent with a Kronecker-product operator
    typedef DefaultMultiVectorProductVector<Scalar> MVPV;
    using Teuchos::dyn_cast;
    using Teuchos::rcp_dynamic_cast;
    RCP<MultiVectorBase<Scalar> > rowStatMultiVec =
        dyn_cast<MVPV>(*rowStatVec).getNonconstMultiVector().assert_not_null();
    const Ordinal numCols = rowStatMultiVec->domain()->dim();
    if (numCols > 0) {
      rcp_dynamic_cast<const RowStatLinearOpBase<Scalar> >(op_.getConstObj())
          ->getRowStat(rowStat, rowStatMultiVec->col(0).ptr());
      for (Ordinal col = 1; col < numCols; ++col) {
        Thyra::copy(*(rowStatMultiVec->col(0)),
                    rowStatMultiVec->col(col).ptr());
      }
    }
  }

  //@}

  /** @name Overridden from ScaledLinearOpBase */
  //@{

  virtual bool supportsScaleLeftImpl() const
  {
    using Teuchos::rcp_dynamic_cast;
    return rcp_dynamic_cast<const ScaledLinearOpBase<Scalar> >(
               op_.getConstObj())
        ->supportsScaleLeft();
  }

  virtual bool supportsScaleRightImpl() const
  {
    using Teuchos::rcp_dynamic_cast;
    return rcp_dynamic_cast<const ScaledLinearOpBase<Scalar> >(
               op_.getConstObj())
        ->supportsScaleRight();
  }

  virtual void scaleLeftImpl(const VectorBase<Scalar> &row_scaling)
  {
    TEUCHOS_ASSERT(this->supportsScaleLeft());

    // Use the first column in the scaling vector to scale the operator
    typedef DefaultMultiVectorProductVector<Scalar> MVPV;
    using Teuchos::dyn_cast;
    using Teuchos::rcp_dynamic_cast;
    RCP<const MultiVectorBase<Scalar> > row_scaling_mv =
        dyn_cast<const MVPV>(row_scaling).getMultiVector().assert_not_null();
    const Ordinal numCols = row_scaling_mv->domain()->dim();
    if (numCols > 0) {
      rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(op_.getNonconstObj())
          ->scaleLeft(*(row_scaling_mv->col(0)));
    }

    // row_scaling.describe(std::cout, Teuchos::VERB_EXTREME);
  }

  virtual void scaleRightImpl(const VectorBase<Scalar> &col_scaling)
  {
    TEUCHOS_ASSERT(this->supportsScaleRight());

    // Use the first column in the scaling vector to scale the operator
    typedef DefaultMultiVectorProductVector<Scalar> MVPV;
    using Teuchos::dyn_cast;
    using Teuchos::rcp_dynamic_cast;
    RCP<const MultiVectorBase<Scalar> > col_scaling_mv =
        dyn_cast<const MVPV>(col_scaling).getMultiVector().assert_not_null();
    const Ordinal numCols = col_scaling_mv->domain()->dim();
    if (numCols > 0) {
      rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(op_.getNonconstObj())
          ->scaleRight(*(col_scaling_mv->col(0)));
    }
  }

  //@}

 private:
  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > CNOP;

  // //////////////////////////////
  // Private data members

  CNOP op_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecRange_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecDomain_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const LinearOpBase<Scalar> > &op,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(op));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecRange));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecDomain));
    TEUCHOS_TEST_FOR_EXCEPT(multiVecRange->numBlocks() !=
                            multiVecDomain->numBlocks());
    if (op->range() != Teuchos::null)
      THYRA_ASSERT_VEC_SPACES(
          "MultiVectorLinearOp<Scalar>::initialize(op,multiVecRange,"
          "multiVecDomain)",
          *op->range(), *multiVecRange->getBlock(0));
    if (op->domain() != Teuchos::null)
      THYRA_ASSERT_VEC_SPACES(
          "MultiVectorLinearOp<Scalar>::initialize(op,multiVecRange,"
          "multiVecDomain)",
          *op->domain(), *multiVecDomain->getBlock(0));
#else
    (void)op;
    (void)multiVecRange;
    (void)multiVecDomain;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorLinearOp
 */
template <class Scalar>
RCP<MultiVectorLinearOp<Scalar> > multiVectorLinearOp()
{
  return Teuchos::rcp(new MultiVectorLinearOp<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorLinearOp
 */
template <class Scalar>
RCP<MultiVectorLinearOp<Scalar> > nonconstMultiVectorLinearOp(
    const RCP<LinearOpBase<Scalar> > &op,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorLinearOp<Scalar> > mvop =
      Teuchos::rcp(new MultiVectorLinearOp<Scalar>());
  mvop->nonconstInitialize(op, multiVecRange, multiVecDomain);
  return mvop;
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorLinearOp
 */
template <class Scalar>
RCP<MultiVectorLinearOp<Scalar> > nonconstMultiVectorLinearOp(
    const RCP<LinearOpBase<Scalar> > &op, const int num_blocks)
{
  RCP<MultiVectorLinearOp<Scalar> > mvop =
      Teuchos::rcp(new MultiVectorLinearOp<Scalar>());
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_domain =
      Thyra::multiVectorProductVectorSpace(op->domain(), num_blocks);
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_range =
      Thyra::multiVectorProductVectorSpace(op->range(), num_blocks);
  mvop->nonconstInitialize(op, mv_range, mv_domain);
  return mvop;
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorLinearOp
 */
template <class Scalar>
RCP<MultiVectorLinearOp<Scalar> > multiVectorLinearOp(
    const RCP<const LinearOpBase<Scalar> > &op,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorLinearOp<Scalar> > mvop =
      Teuchos::rcp(new MultiVectorLinearOp<Scalar>());
  mvop->initialize(op, multiVecRange, multiVecDomain);
  return mvop;
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorLinearOp
 */
template <class Scalar>
RCP<MultiVectorLinearOp<Scalar> > multiVectorLinearOp(
    const RCP<const LinearOpBase<Scalar> > &op, const int num_blocks)
{
  RCP<MultiVectorLinearOp<Scalar> > mvop =
      Teuchos::rcp(new MultiVectorLinearOp<Scalar>());
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_domain =
      Thyra::multiVectorProductVectorSpace(op->domain(), num_blocks);
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_range =
      Thyra::multiVectorProductVectorSpace(op->range(), num_blocks);
  mvop->initialize(op, mv_range, mv_domain);
  return mvop;
}

}  // end namespace Thyra

#endif
