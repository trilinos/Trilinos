// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockDiagonalInverseOp.hpp"

#include "Teuchos_Utils.hpp"

namespace Teko {

using Teuchos::RCP;

BlockDiagonalInverseOp::BlockDiagonalInverseOp(BlockedLinearOp& A,
                                               const std::vector<LinearOp>& invDiag)
    : invDiag_(invDiag) {
  // sanity check
  int blocks = blockRowCount(A);
  TEUCHOS_ASSERT(blocks > 0);
  TEUCHOS_ASSERT(blocks == blockColCount(A));
  TEUCHOS_ASSERT(blocks == (int)invDiag_.size());

  // create the range and product space
  ///////////////////////////////////////////////////

  // just flip flop them!
  productRange_  = A->productDomain();
  productDomain_ = A->productRange();
}

void BlockDiagonalInverseOp::implicitApply(const BlockedMultiVector& src, BlockedMultiVector& dst,
                                           const double alpha, const double beta) const {
  // call the no tranpose version
  implicitApply(Thyra::NOTRANS, src, dst, alpha, beta);
}

void BlockDiagonalInverseOp::implicitApply(const Thyra::EOpTransp M_trans,
                                           const BlockedMultiVector& src, BlockedMultiVector& dst,
                                           const double alpha, const double beta) const {
  int blocks = blockCount(src);

  TEUCHOS_ASSERT(blocks == (int)invDiag_.size());

  if (!allocated) {
    srcScrap_ = deepcopy(src);
    dstScrap_ = deepcopy(dst);
    allocated = true;
  }

  // build a scrap vector for storing work
  Thyra::assign<double>(srcScrap_.ptr(), *src);
  BlockedMultiVector dstCopy;
  if (beta != 0.0) {
    Thyra::assign<double>(dstScrap_.ptr(), *dst);
    dstCopy = dstScrap_;
  } else
    dstCopy = dst;  // shallow copy

  // extract the blocks components from
  // the source and destination vectors
  std::vector<MultiVector> dstVec;
  std::vector<MultiVector> scrapVec;
  for (int b = 0; b < blocks; b++) {
    dstVec.push_back(getBlock(b, dstCopy));
    scrapVec.push_back(getBlock(b, srcScrap_));
  }

  if (M_trans == Thyra::NOTRANS) {
    for (int b = 0; b < blocks; b++) {
      applyOp(invDiag_[b], scrapVec[b], dstVec[b]);
    }
  } else if (M_trans == Thyra::TRANS || M_trans == Thyra::CONJTRANS) {
    for (int b = 0; b < blocks; b++) {
      applyTransposeOp(invDiag_[b], scrapVec[b], dstVec[b]);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  // scale result by alpha
  if (beta != 0)
    update(alpha, dstCopy, beta, dst);  // dst = alpha * dstCopy + beta * dst
  else if (alpha != 1.0)
    scale(alpha, dst);  // dst = alpha * dst
}

void BlockDiagonalInverseOp::describe(Teuchos::FancyOStream& out_arg,
                                      const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::OSTab;

  RCP<Teuchos::FancyOStream> out = rcp(&out_arg, false);
  OSTab tab(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW: *out << this->description() << std::endl; break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      *out << Teuchos::Describable::description() << "{"
           << "rangeDim=" << this->range()->dim() << ",domainDim=" << this->domain()->dim()
           << ",rows=" << invDiag_.size() << ",cols=" << invDiag_.size() << "}\n";
      {
        OSTab tab2(out);
        *out << "[invDiag Operators]:\n";
        tab.incrTab();
        for (auto i = 0U; i < invDiag_.size(); i++) {
          *out << "[invD(" << i << ")] = ";
          *out << Teuchos::describe(*invDiag_[i], verbLevel);
        }
      }
      break;
    }
    default: TEUCHOS_TEST_FOR_EXCEPT(true);  // Should never get here!
  }
}

}  // end namespace Teko
