// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_BLOCKED_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultBlockedLinearOp_decl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"

namespace Thyra {

// Constructors

template <class Scalar>
DefaultBlockedLinearOp<Scalar>::DefaultBlockedLinearOp()
  : numRowBlocks_(0)
  , numColBlocks_(0)
  , blockFillIsActive_(false) {}

// Overridden from PhysicallyBlockedLinearOpBase

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill() {
  assertBlockFillIsActive(false);
  uninitialize();
  resetStorage(0, 0);
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
    const int numRowBlocks, const int numColBlocks) {
  assertBlockFillIsActive(false);
  uninitialize();
  resetStorage(numRowBlocks, numColBlocks);
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
    const RCP<const ProductVectorSpaceBase<Scalar> > &new_productRange, const RCP<const ProductVectorSpaceBase<Scalar> > &new_productDomain) {
  using Teuchos::rcp_dynamic_cast;
  assertBlockFillIsActive(false);
  uninitialize();
  productRange_  = new_productRange.assert_not_null();
  productDomain_ = new_productDomain.assert_not_null();
  defaultProductRange_ =
      rcp_dynamic_cast<const DefaultProductVectorSpace<Scalar> >(productRange_);
  defaultProductDomain_ =
      rcp_dynamic_cast<const DefaultProductVectorSpace<Scalar> >(productDomain_);
  // Note: the above spaces must be set before calling the next function!
  resetStorage(productRange_->numBlocks(), productDomain_->numBlocks());
}

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockFillIsActive() const {
  return blockFillIsActive_;
}

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::acceptsBlock(
    const int i, const int j) const {
  assertBlockFillIsActive(true);
  assertBlockRowCol(i, j);
  return true;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::setNonconstBlock(
    const int i, const int j, const RCP<LinearOpBase<Scalar> > &block) {
  setBlockImpl(i, j, block);
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlock(
    const int i, const int j, const RCP<const LinearOpBase<Scalar> > &block) {
  setBlockImpl(i, j, block);
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::endBlockFill() {
  using Teuchos::as;

  assertBlockFillIsActive(true);

  // 2009/05/06: rabartl: ToDo: When doing a flexible block fill
  // (Ops_stack_.size() > 0), we need to assert that all of the block rows and
  // columns have been filled in.  I don't think we do that here.

  // Get the number of block rows and columns
  if (nonnull(productRange_)) {
    numRowBlocks_ = productRange_->numBlocks();
    numColBlocks_ = productDomain_->numBlocks();
  } else {
    numRowBlocks_ = rangeBlocks_.size();
    numColBlocks_ = domainBlocks_.size();
    // NOTE: Above, whether doing a flexible fill or not, all of the blocks
    // must be set in order to have a valid filled operator so this
    // calculation should be correct.
  }

  // Assert that all of the block rows and columns have at least one entry if
  // the spaces were not given up front.
#ifdef TEUCHOS_DEBUG
  if (is_null(productRange_)) {
    for (int i = 0; i < numRowBlocks_; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          !rangeBlocks_[i].get(), std::logic_error,
          "DefaultBlockedLinearOp<Scalar>::endBlockFill():"
          " Error, no linear operator block for the i="
              << i << " block row was added"
                      " and we can not complete the block fill!");
    }
    for (int j = 0; j < numColBlocks_; ++j) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          !domainBlocks_[j].get(), std::logic_error,
          "DefaultBlockedLinearOp<Scalar>::endBlockFill():"
          " Error, no linear operator block for the j="
              << j << " block column was added"
                      " and we can not complete the block fill!");
    }
  }
#endif

  // Insert the block LOB objects if doing a flexible fill.
  if (Ops_stack_.size()) {
    Ops_.resize(numRowBlocks_ * numColBlocks_);
    for (int k = 0; k < as<int>(Ops_stack_.size()); ++k) {
      const BlockEntry<Scalar> &block_i_j             = Ops_stack_[k];
      Ops_[numRowBlocks_ * block_i_j.j + block_i_j.i] = block_i_j.block;
    }
    Ops_stack_.resize(0);
  }

  // Set the product range and domain spaces if not already set
  if (is_null(productRange_)) {
    adjustBlockSpaces();
    defaultProductRange_  = productVectorSpace<Scalar>(rangeBlocks_());
    defaultProductDomain_ = productVectorSpace<Scalar>(domainBlocks_());
    productRange_         = defaultProductRange_;
    productDomain_        = defaultProductDomain_;
  }

  rangeBlocks_.resize(0);
  domainBlocks_.resize(0);

  blockFillIsActive_ = false;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::uninitialize() {
  productRange_  = Teuchos::null;
  productDomain_ = Teuchos::null;
  numRowBlocks_  = 0;
  numColBlocks_  = 0;
  Ops_.resize(0);
  Ops_stack_.resize(0);
  rangeBlocks_.resize(0);
  domainBlocks_.resize(0);
  blockFillIsActive_ = false;
}

// Overridden from BlockedLinearOpBase

template <class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productRange() const {
  return productRange_;
}

template <class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productDomain() const {
  return productDomain_;
}

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockExists(
    const int i, const int j) const {
  assertBlockFillIsActive(false);
  assertBlockRowCol(i, j);
  return true;
}

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockIsConst(
    const int i, const int j) const {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i, j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i, j);
  return Ops_[numRowBlocks_ * j + i].isConst();
}

template <class Scalar>
RCP<LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getNonconstBlock(const int i, const int j) {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i, j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i, j);
  return Ops_[numRowBlocks_ * j + i].getNonconstObj();
}

template <class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getBlock(const int i, const int j) const {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i, j));
#endif
  assertBlockFillIsActive(false);
  assertBlockRowCol(i, j);
  return Ops_[numRowBlocks_ * j + i];
}

// Overridden from LinearOpBase

template <class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::range() const {
  return productRange_;
}

template <class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::domain() const {
  return productDomain_;
}

template <class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::clone() const {
  return Teuchos::null;  // ToDo: Implement this when needed!
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string DefaultBlockedLinearOp<Scalar>::description() const {
  assertBlockFillIsActive(false);
  std::ostringstream oss;
  oss
      << Teuchos::Describable::description() << "{"
      << "numRowBlocks=" << numRowBlocks_
      << ",numColBlocks=" << numColBlocks_
      << "}";
  return oss.str();
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::describe(
    Teuchos::FancyOStream &out_arg, const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::rcpFromRef;
  assertBlockFillIsActive(false);
  RCP<FancyOStream> out = rcpFromRef(out_arg);
  OSTab tab1(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      *out
          << Teuchos::Describable::description() << "{"
          << "rangeDim=" << this->range()->dim()
          << ",domainDim=" << this->domain()->dim()
          << ",numRowBlocks=" << numRowBlocks_
          << ",numColBlocks=" << numColBlocks_
          << "}\n";
      OSTab tab2(out);
      *out
          << "Constituent LinearOpBase objects for M = [ Op[0,0] ..."
          << " ; ... ; ... Op[numRowBlocks-1,numColBlocks-1] ]:\n";
      tab2.incrTab();
      for (int i = 0; i < numRowBlocks_; ++i) {
        for (int j = 0; j < numColBlocks_; ++j) {
          *out << "Op[" << i << "," << j << "] = ";
          RCP<const LinearOpBase<Scalar> >
              block_i_j = getBlock(i, j);
          if (block_i_j.get())
            *out << Teuchos::describe(*getBlock(i, j), verbLevel);
          else
            *out << "NULL\n";
        }
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);  // Should never get here!
  }
}

// protected

// Overridden from LinearOpBase

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const {
  bool supported = true;
  for (int i = 0; i < numRowBlocks_; ++i) {
    for (int j = 0; j < numColBlocks_; ++j) {
      RCP<const LinearOpBase<Scalar> >
          block_i_j = getBlock(i, j);
      if (block_i_j.get() && !Thyra::opSupported(*block_i_j, M_trans))
        supported = false;
    }
  }
  return supported;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X_in,
    const Ptr<MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta) const {
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef RCP<MultiVectorBase<Scalar> > MultiVectorPtr;
  typedef RCP<const MultiVectorBase<Scalar> > ConstMultiVectorPtr;
  typedef RCP<const LinearOpBase<Scalar> > ConstLinearOpPtr;

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
      "DefaultBlockedLinearOp<Scalar>::apply(...)", *this, M_trans, X_in, &*Y_inout);
#endif  // TEUCHOS_DEBUG

  const bool
      struct_transp = (real_trans(M_trans) != NOTRANS);  // Structural transpose?
  const int
      opNumRowBlocks = (!struct_transp ? numRowBlocks_ : numColBlocks_),
      opNumColBlocks = (!struct_transp ? numColBlocks_ : numRowBlocks_);

  //
  // Y = alpha * op(M) * X + beta*Y
  //
  // =>
  //
  // Y[i] = beta+Y[i] + sum(alpha*op(Op)[i,j]*X[j],j=0...opNumColBlocks-1)
  //
  // , for i=0...opNumRowBlocks-1
  //

  const RCP<const DefaultProductVectorSpace<Scalar> >
      defaultProductRange_op  = (real_trans(M_trans) == NOTRANS
                                     ? defaultProductRange_
                                     : defaultProductDomain_),
      defaultProductDomain_op = (real_trans(M_trans) == NOTRANS
                                     ? defaultProductDomain_
                                     : defaultProductRange_);

  const RCP<const ProductMultiVectorBase<Scalar> >
      X = castOrCreateSingleBlockProductMultiVector<Scalar>(
          defaultProductDomain_op, rcpFromRef(X_in));
  const RCP<ProductMultiVectorBase<Scalar> >
      Y = nonconstCastOrCreateSingleBlockProductMultiVector<Scalar>(
          defaultProductRange_op, rcpFromPtr(Y_inout));

  for (int i = 0; i < opNumRowBlocks; ++i) {
    MultiVectorPtr Y_i = Y->getNonconstMultiVectorBlock(i);
    for (int j = 0; j < opNumColBlocks; ++j) {
      ConstLinearOpPtr
          Op_i_j = (!struct_transp ? getBlock(i, j) : getBlock(j, i));
      ConstMultiVectorPtr
          X_j = X->getMultiVectorBlock(j);
      if (j == 0) {
        if (nonnull(Op_i_j))
          Thyra::apply(*Op_i_j, M_trans, *X_j, Y_i.ptr(), alpha, beta);
        else
          scale(beta, Y_i.ptr());
      } else {
        if (nonnull(Op_i_j))
          Thyra::apply(*Op_i_j, M_trans, *X_j, Y_i.ptr(), alpha, ST::one());
      }
    }
  }
}

// Overridden from RowStatLinearOpBase

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::rowStatIsSupportedImpl(
    const RowStatLinearOpBaseUtils::ERowStat row_stat) const {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  // typedef Teuchos::ScalarTraits<Scalar> ST; // unused
  typedef RowStatLinearOpBase<Scalar> RowStatOp;
  typedef RCP<const LinearOpBase<Scalar> > ConstLinearOpPtr;
  typedef const LinearOpBase<Scalar> ConstLinearOp;

  // Figure out what needs to be check for each sub block to compute
  // the require row stat
  RowStatLinearOpBaseUtils::ERowStat
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
  switch (row_stat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM;
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  // check each sub block for compatibility (null means zero,
  // automatically compatible)
  for (int i = 0; i < numRowBlocks_; ++i) {
    for (int j = 0; j < numColBlocks_; ++j) {
      ConstLinearOpPtr Op_i_j = getBlock(i, j);

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        ConstLinearOp *Orig_i_j = 0;
        Scalar scalar           = 1.0;
        EOpTransp transp        = NOTRANS;
        Thyra::unwrap(*Op_i_j, &scalar, &transp, &Orig_i_j);

        const RowStatOp &row_stat_op = Teuchos::dyn_cast<const RowStatOp>(*Orig_i_j);

        // sub block must also support the required row stat operation
        RowStatLinearOpBaseUtils::ERowStat stat = subblk_stat;
        if (transp == NOTRANS || transp == CONJ)
          stat = subblk_stat;
        else if (transp == TRANS || transp == CONJTRANS)
          stat = subblk_trans_stat;

        if (!row_stat_op.rowStatIsSupported(stat))
          return false;
      }
      // else: sub block is null, "0" is good enough
    }
  }

  return true;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat row_stat,
    const Teuchos::Ptr<VectorBase<Scalar> > &rowStatVec) const {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromPtr;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef RowStatLinearOpBase<Scalar> RowStatOp;
  typedef RCP<const LinearOpBase<Scalar> > ConstLinearOpPtr;
  typedef const LinearOpBase<Scalar> ConstLinearOp;
  typedef RCP<VectorBase<Scalar> > VectorPtr;

  // Figure out what needs to be check for each sub block to compute
  // the require row stat
  RowStatLinearOpBaseUtils::ERowStat
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
  switch (row_stat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM;
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      subblk_stat       = RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
      subblk_trans_stat = RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  const RCP<ProductVectorBase<Scalar> >
      Y = rcp_dynamic_cast<ProductVectorBase<Scalar> >(rcpFromPtr(rowStatVec));

  // using sub block row stat capability interrogate each for
  // its constribution
  for (int i = 0; i < numRowBlocks_; ++i) {
    VectorPtr blk_vec = Y->getNonconstVectorBlock(i);
    put_scalar(ST::zero(), blk_vec.ptr());  // clear vector just in case

    for (int j = 0; j < numColBlocks_; ++j) {
      ConstLinearOpPtr Op_i_j = getBlock(i, j);

      VectorPtr tmp_vec = createMember(Op_i_j->range());

      put_scalar(ST::zero(), tmp_vec.ptr());  // clear vector just in case

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        ConstLinearOp *Orig_i_j = 0;
        Scalar scalar           = 1.0;
        EOpTransp transp        = NOTRANS;
        Thyra::unwrap(*Op_i_j, &scalar, &transp, &Orig_i_j);

        const RowStatOp &row_stat_op = Teuchos::dyn_cast<const RowStatOp>(*Orig_i_j);

        // sub block must also support the required row stat operation
        if (transp == NOTRANS || transp == CONJ)
          row_stat_op.getRowStat(subblk_stat, tmp_vec.ptr());
        else if (transp == TRANS || transp == CONJTRANS)
          row_stat_op.getRowStat(subblk_trans_stat, tmp_vec.ptr());
        else {
          TEUCHOS_ASSERT(false);
        }

        // some the temporary into the block vector
        Vp_V(blk_vec.ptr(), *tmp_vec);
      }
    }
  }

  // do post processing as needed (take the inverse currently
  switch (row_stat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
      reciprocal(*rowStatVec, rowStatVec.ptr());
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}

// Overridden from ScaledLinearOpBase

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::supportsScaleLeftImpl() const {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef RCP<const LinearOpBase<Scalar> > LinearOpPtr;
  typedef const LinearOpBase<Scalar> LinearOp;
  typedef const ScaledLinearOpBase<Scalar> ScaledOp;

  bool supported = true;
  for (int i = 0; i < numRowBlocks_; ++i) {
    for (int j = 0; j < numColBlocks_; ++j) {
      LinearOpPtr Op_i_j = getBlock(i, j);

      EOpTransp transp = NOTRANS;

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        LinearOp *Orig_i_j = 0;
        Scalar scalar      = 1.0;
        Thyra::unwrap(*Op_i_j, &scalar, &transp, &Orig_i_j);

        ScaledOp &scaled_op = Teuchos::dyn_cast<ScaledOp>(*Orig_i_j);

        if (transp == NOTRANS || transp == CONJ)
          supported &= scaled_op.supportsScaleLeft();
        else if (transp == TRANS || transp == CONJTRANS)
          supported &= scaled_op.supportsScaleRight();
      }
    }
  }

  return supported;
}

template <class Scalar>
bool DefaultBlockedLinearOp<Scalar>::supportsScaleRightImpl() const {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef RCP<const LinearOpBase<Scalar> > LinearOpPtr;
  typedef const LinearOpBase<Scalar> LinearOp;
  typedef const ScaledLinearOpBase<Scalar> ScaledOp;

  bool supported = true;
  for (int i = 0; i < numRowBlocks_; ++i) {
    for (int j = 0; j < numColBlocks_; ++j) {
      LinearOpPtr Op_i_j = getBlock(i, j);

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        LinearOp *Orig_i_j = 0;
        Scalar scalar      = 1.0;
        EOpTransp transp   = NOTRANS;
        Thyra::unwrap(*Op_i_j, &scalar, &transp, &Orig_i_j);

        ScaledOp &scaled_op = Teuchos::dyn_cast<ScaledOp>(*Orig_i_j);

        if (transp == NOTRANS || transp == CONJ)
          supported &= scaled_op.supportsScaleRight();
        else if (transp == TRANS || transp == CONJTRANS)
          supported &= scaled_op.supportsScaleLeft();
      }
    }
  }

  return supported;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::scaleLeftImpl(
    const VectorBase<Scalar> &row_scaling) {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef ScaledLinearOpBase<Scalar> ScaledOp;
  typedef RCP<LinearOpBase<Scalar> > LinearOpPtr;
  typedef RCP<const VectorBase<Scalar> > VectorPtr;

  const ProductVectorBase<Scalar> &
      Y = dyn_cast<const ProductVectorBase<Scalar> >(row_scaling);

  // using sub block row stat capability interrogate each for
  // its constribution
  for (int i = 0; i < numRowBlocks_; ++i) {
    VectorPtr blk_vec = Y.getVectorBlock(i);

    for (int j = 0; j < numColBlocks_; ++j) {
      LinearOpPtr Op_i_j = getNonconstBlock(i, j);

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        LinearOpPtr Orig_i_j;
        EOpTransp transp = NOTRANS;

        {
          RCP<ScaledAdjointLinearOpBase<Scalar> > saOp = rcp_dynamic_cast<ScaledAdjointLinearOpBase<Scalar> >(Op_i_j);
          if (nonnull(saOp)) {
            transp   = saOp->overallTransp();
            Orig_i_j = saOp->getNonconstOrigOp();
          } else
            Orig_i_j = Op_i_j;  // stick with the original
        }

        RCP<ScaledOp> scaled_op = rcp_dynamic_cast<ScaledOp>(Orig_i_j);

        // sub block must support a row stat operation
        TEUCHOS_ASSERT(scaled_op != Teuchos::null);  // guranteed by compatibility check

        // sub block must also support the required row stat operation
        if (transp == NOTRANS || transp == CONJ)
          scaled_op->scaleLeft(*blk_vec);
        else if (transp == TRANS || transp == CONJTRANS)
          scaled_op->scaleRight(*blk_vec);
      }
    }
  }
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::scaleRightImpl(
    const VectorBase<Scalar> &col_scaling) {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef ScaledLinearOpBase<Scalar> ScaledOp;
  typedef RCP<LinearOpBase<Scalar> > LinearOpPtr;
  typedef RCP<const VectorBase<Scalar> > VectorPtr;

  const ProductVectorBase<Scalar> &
      Y = dyn_cast<const ProductVectorBase<Scalar> >(col_scaling);

  // using sub block row stat capability interrogate each for
  // its constribution
  for (int j = 0; j < numColBlocks_; ++j) {
    VectorPtr blk_vec = Y.getVectorBlock(j);

    for (int i = 0; i < numRowBlocks_; ++i) {
      LinearOpPtr Op_i_j = getNonconstBlock(i, j);

      if (nonnull(Op_i_j)) {
        // pull out adjoint, and scaling
        LinearOpPtr Orig_i_j;
        EOpTransp transp = NOTRANS;

        {
          RCP<ScaledAdjointLinearOpBase<Scalar> > saOp = rcp_dynamic_cast<ScaledAdjointLinearOpBase<Scalar> >(Op_i_j);
          if (nonnull(saOp)) {
            transp   = saOp->overallTransp();
            Orig_i_j = saOp->getNonconstOrigOp();
          } else
            Orig_i_j = Op_i_j;  // stick with the original
        }

        RCP<ScaledOp> scaled_op = rcp_dynamic_cast<ScaledOp>(Orig_i_j);

        // sub block must support a row stat operation
        TEUCHOS_ASSERT(scaled_op != Teuchos::null);  // guranteed by compatibility check

        // sub block must also support the required row stat operation
        if (transp == NOTRANS || transp == CONJ)
          scaled_op->scaleRight(*blk_vec);
        else if (transp == TRANS || transp == CONJTRANS)
          scaled_op->scaleLeft(*blk_vec);
      }
    }
  }
}

// private

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::resetStorage(
    const int numRowBlocks, const int numColBlocks) {
  numRowBlocks_ = numRowBlocks;
  numColBlocks_ = numColBlocks;
  Ops_.resize(numRowBlocks_ * numColBlocks_);
  if (is_null(productRange_)) {
    rangeBlocks_.resize(numRowBlocks);
    domainBlocks_.resize(numColBlocks);
  }
  blockFillIsActive_ = true;
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::assertBlockFillIsActive(
    bool wantedValue) const {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!(blockFillIsActive_ == wantedValue));
#else
  (void)wantedValue;
#endif
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::assertBlockRowCol(
    const int i, const int j) const {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(0 <= i), std::logic_error, "Error, i=" << i << " is invalid!");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(0 <= j), std::logic_error, "Error, j=" << j << " is invalid!");
  // Only validate upper range if the number of row and column blocks is
  // fixed!
  if (Ops_.size()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !(0 <= i && i < numRowBlocks_), std::logic_error, "Error, i=" << i << " does not fall in the range [0," << numRowBlocks_ - 1 << "]!");
    TEUCHOS_TEST_FOR_EXCEPTION(
        !(0 <= j && j < numColBlocks_), std::logic_error, "Error, j=" << j << " does not fall in the range [0," << numColBlocks_ - 1 << "]!");
  }
#else
  (void)i;
  (void)j;
#endif
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlockSpaces(
    const int i, const int j, const LinearOpBase<Scalar> &block) {
  using Teuchos::toString;
  assertBlockFillIsActive(true);
  assertBlockRowCol(i, j);

  // Validate that if the vector space block is already set that it is
  // compatible with the block that is being set.
  if (i < numRowBlocks_ && j < numColBlocks_) {
#ifdef TEUCHOS_DEBUG
    RCP<const VectorSpaceBase<Scalar> >
        rangeBlock  = (productRange_.get()
                           ? productRange_->getBlock(i)
                           : rangeBlocks_[i]),
        domainBlock = (productDomain_.get()
                           ? productDomain_->getBlock(j)
                           : domainBlocks_[j]);
    if (rangeBlock.get()) {
      THYRA_ASSERT_VEC_SPACES_NAMES(
          "DefaultBlockedLinearOp<Scalar>::setBlockSpaces(i,j,block):\n\n"
          "Adding block: " +
              block.description(),
          *rangeBlock, ("(*productRange->getBlock(" + toString(i) + "))"),
          *block.range(), ("(*block[" + toString(i) + "," + toString(j) + "].range())"));
    }
    if (domainBlock.get()) {
      THYRA_ASSERT_VEC_SPACES_NAMES(
          "DefaultBlockedLinearOp<Scalar>::setBlockSpaces(i,j,block):\n\n"
          "Adding block: " +
              block.description(),
          *domainBlock, ("(*productDomain->getBlock(" + toString(j) + "))"),
          *block.domain(), ("(*block[" + toString(i) + "," + toString(j) + "].domain())"));
    }
#endif  // TEUCHOS_DEBUG
  }

  // Add spaces missing range and domain space blocks if we are doing a
  // flexible fill (otherwise these loops will not be executed)
  for (int k = numRowBlocks_; k <= i; ++k)
    rangeBlocks_.push_back(Teuchos::null);
  for (int k = numColBlocks_; k <= j; ++k)
    domainBlocks_.push_back(Teuchos::null);

  // Set the incoming range and domain blocks if not already set
  if (!productRange_.get()) {
    if (!rangeBlocks_[i].get())
      rangeBlocks_[i] = block.range().assert_not_null();
    if (!domainBlocks_[j].get()) {
      domainBlocks_[j] = block.domain().assert_not_null();
    }
  }

  // Update the current number of row and columns blocks if doing a flexible
  // fill.
  if (!Ops_.size()) {
    numRowBlocks_ = rangeBlocks_.size();
    numColBlocks_ = domainBlocks_.size();
  }
}

template <class Scalar>
template <class LinearOpType>
void DefaultBlockedLinearOp<Scalar>::setBlockImpl(
    const int i, const int j,
    const RCP<LinearOpType> &block) {
  setBlockSpaces(i, j, *block);
  if (Ops_.size()) {
    // We are doing a fill with a fixed number of row and column blocks so we
    // can just set this.
    Ops_[numRowBlocks_ * j + i] = block;
  } else {
    // We are doing a flexible fill so add the block to the stack of blocks or
    // replace a block that already exists.
    bool foundBlock = false;
    for (unsigned int k = 0; k < Ops_stack_.size(); ++k) {
      BlockEntry<Scalar> &block_i_j = Ops_stack_[k];
      if (block_i_j.i == i && block_i_j.j == j) {
        block_i_j.block = block;
        foundBlock      = true;
        break;
      }
    }
    if (!foundBlock)
      Ops_stack_.push_back(BlockEntry<Scalar>(i, j, block));
  }
}

template <class Scalar>
void DefaultBlockedLinearOp<Scalar>::adjustBlockSpaces() {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INEQUALITY(Ops_.size(), !=, 0);
#endif

  //
  // Loop through the rows and columns looking for rows with mixed
  // single-space range and/or domain spaces on operators and set the single
  // spaces as the block space if it exists.
  //
  // NOTE: Once we get here, we can safely assume that all of the operators
  // are compatible w.r.t. their spaces so if there are rows and/or columns
  // with mixed product and single vector spaces that we can just pick the
  // single vector space for the whole row and/or column.
  //

  // Adjust blocks in the range space
  for (int i = 0; i < numRowBlocks_; ++i) {
    for (int j = 0; j < numColBlocks_; ++j) {
      const RCP<const LinearOpBase<Scalar> >
          op_i_j = Ops_[numRowBlocks_ * j + i];
      if (is_null(op_i_j))
        continue;
      const RCP<const VectorSpaceBase<Scalar> > range_i_j = op_i_j->range();
      if (is_null(productVectorSpaceBase<Scalar>(range_i_j, false))) {
        rangeBlocks_[i] = range_i_j;
        break;
      }
    }
  }

  // Adjust blocks in the domain space
  for (int j = 0; j < numColBlocks_; ++j) {
    for (int i = 0; i < numRowBlocks_; ++i) {
      const RCP<const LinearOpBase<Scalar> >
          op_i_j = Ops_[numRowBlocks_ * j + i];
      if (is_null(op_i_j))
        continue;
      const RCP<const VectorSpaceBase<Scalar> >
          domain_i_j = op_i_j->domain();
      if (is_null(productVectorSpaceBase<Scalar>(domain_i_j, false))) {
        domainBlocks_[j] = domain_i_j;
        break;
      }
    }
  }
}

}  // namespace Thyra

template <class Scalar>
Teuchos::RCP<Thyra::DefaultBlockedLinearOp<Scalar> >
Thyra::defaultBlockedLinearOp() {
  return Teuchos::rcp(new DefaultBlockedLinearOp<Scalar>());
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::block1x1(
    const RCP<const LinearOpBase<Scalar> > &A00,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(1, 1);
  M->setBlock(0, 0, A00);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::block1x2(
    const RCP<const LinearOpBase<Scalar> > &A00,
    const RCP<const LinearOpBase<Scalar> > &A01,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(1, 2);
  if (A00.get()) M->setBlock(0, 0, A00);
  if (A01.get()) M->setBlock(0, 1, A01);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::block2x1(
    const RCP<const LinearOpBase<Scalar> > &A00,
    const RCP<const LinearOpBase<Scalar> > &A10,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(2, 1);
  if (A00.get()) M->setBlock(0, 0, A00);
  if (A10.get()) M->setBlock(1, 0, A10);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::block2x2(
    const RCP<const LinearOpBase<Scalar> > &A00,
    const RCP<const LinearOpBase<Scalar> > &A01,
    const RCP<const LinearOpBase<Scalar> > &A10,
    const RCP<const LinearOpBase<Scalar> > &A11,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(2, 2);
  if (A00.get()) M->setBlock(0, 0, A00);
  if (A01.get()) M->setBlock(0, 1, A01);
  if (A10.get()) M->setBlock(1, 0, A10);
  if (A11.get()) M->setBlock(1, 1, A11);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock1x1(
    const RCP<LinearOpBase<Scalar> > &A00,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(1, 1);
  M->setNonconstBlock(0, 0, A00);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock1x2(
    const RCP<LinearOpBase<Scalar> > &A00,
    const RCP<LinearOpBase<Scalar> > &A01,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(1, 2);
  if (A00.get()) M->setNonconstBlock(0, 0, A00);
  if (A01.get()) M->setNonconstBlock(0, 1, A01);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock2x1(
    const RCP<LinearOpBase<Scalar> > &A00,
    const RCP<LinearOpBase<Scalar> > &A10,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(2, 1);
  if (A00.get()) M->setNonconstBlock(0, 0, A00);
  if (A10.get()) M->setNonconstBlock(1, 0, A10);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstBlock2x2(
    const RCP<LinearOpBase<Scalar> > &A00,
    const RCP<LinearOpBase<Scalar> > &A01,
    const RCP<LinearOpBase<Scalar> > &A10,
    const RCP<LinearOpBase<Scalar> > &A11,
    const std::string &label) {
  RCP<PhysicallyBlockedLinearOpBase<Scalar> >
      M = defaultBlockedLinearOp<Scalar>();
  M->beginBlockFill(2, 2);
  if (A00.get()) M->setNonconstBlock(0, 0, A00);
  if (A01.get()) M->setNonconstBlock(0, 1, A01);
  if (A10.get()) M->setNonconstBlock(1, 0, A10);
  if (A11.get()) M->setNonconstBlock(1, 1, A11);
  M->endBlockFill();
  if (label.length())
    M->setObjectLabel(label);
  return M;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//

#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_INSTANT(SCALAR) \
                                                        \
  template class DefaultBlockedLinearOp<SCALAR>;        \
                                                        \
  template RCP<DefaultBlockedLinearOp<SCALAR> >         \
  defaultBlockedLinearOp<SCALAR>();                     \
                                                        \
  template RCP<const LinearOpBase<SCALAR> >             \
  block1x1(                                             \
      const RCP<const LinearOpBase<SCALAR> > &A00,      \
      const std::string &label);                        \
                                                        \
  template RCP<const LinearOpBase<SCALAR> >             \
  block1x2(                                             \
      const RCP<const LinearOpBase<SCALAR> > &A00,      \
      const RCP<const LinearOpBase<SCALAR> > &A01,      \
      const std::string &label);                        \
                                                        \
  template RCP<const LinearOpBase<SCALAR> >             \
  block2x1(                                             \
      const RCP<const LinearOpBase<SCALAR> > &A00,      \
      const RCP<const LinearOpBase<SCALAR> > &A10,      \
      const std::string &label);                        \
                                                        \
  template RCP<const LinearOpBase<SCALAR> >             \
  block2x2(                                             \
      const RCP<const LinearOpBase<SCALAR> > &A00,      \
      const RCP<const LinearOpBase<SCALAR> > &A01,      \
      const RCP<const LinearOpBase<SCALAR> > &A10,      \
      const RCP<const LinearOpBase<SCALAR> > &A11,      \
      const std::string &label);                        \
                                                        \
  template RCP<LinearOpBase<SCALAR> >                   \
  nonconstBlock1x1(                                     \
      const RCP<LinearOpBase<SCALAR> > &A00,            \
      const std::string &label);                        \
                                                        \
  template RCP<LinearOpBase<SCALAR> >                   \
  nonconstBlock1x2(                                     \
      const RCP<LinearOpBase<SCALAR> > &A00,            \
      const RCP<LinearOpBase<SCALAR> > &A01,            \
      const std::string &label);                        \
                                                        \
  template RCP<LinearOpBase<SCALAR> >                   \
  nonconstBlock2x1(                                     \
      const RCP<LinearOpBase<SCALAR> > &A00,            \
      const RCP<LinearOpBase<SCALAR> > &A10,            \
      const std::string &label);                        \
                                                        \
  template RCP<LinearOpBase<SCALAR> >                   \
  nonconstBlock2x2(                                     \
      const RCP<LinearOpBase<SCALAR> > &A00,            \
      const RCP<LinearOpBase<SCALAR> > &A01,            \
      const RCP<LinearOpBase<SCALAR> > &A10,            \
      const RCP<LinearOpBase<SCALAR> > &A11,            \
      const std::string &label);

#endif  // THYRA_DEFAULT_BLOCKED_LINEAR_OP_DEF_HPP
