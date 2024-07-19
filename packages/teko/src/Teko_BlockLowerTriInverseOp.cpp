// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockLowerTriInverseOp.hpp"

#include "Teuchos_Utils.hpp"

namespace Teko {

using Teuchos::RCP;

BlockLowerTriInverseOp::BlockLowerTriInverseOp(BlockedLinearOp& L,
                                               const std::vector<LinearOp>& invDiag)
    : L_(L) {
  invDiag_ = invDiag;

  // sanity check
  int blocks = blockRowCount(L_);
  TEUCHOS_ASSERT(blocks > 0);
  TEUCHOS_ASSERT(blocks == blockColCount(L_));
  TEUCHOS_ASSERT(blocks == (int)invDiag_.size());

  // create the range and product space
  ///////////////////////////////////////////////////

  // just flip flop them!
  productRange_  = L_->productDomain();
  productDomain_ = L_->productRange();
}

/** @brief Perform a matrix vector multiply with this operator.
 *
 * The <code>apply</code> function takes one vector as input
 * and applies the inverse \f$ LDU \f$ decomposition. The result
 * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
 * \f$ y = \alpha M x + \beta y \f$ (ignoring conjugation!).
 *
 * @param[in]     x
 * @param[in,out] y
 * @param[in]     alpha (default=1)
 * @param[in]     beta  (default=0)
 */
void BlockLowerTriInverseOp::implicitApply(const BlockedMultiVector& src, BlockedMultiVector& dst,
                                           const double alpha, const double beta) const {
  int blocks = blockCount(src);

  TEUCHOS_ASSERT(blocks == blockRowCount(L_));
  TEUCHOS_ASSERT(blocks == blockCount(dst));

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

  // extract the blocks componets from
  // the source and destination vectors
  std::vector<MultiVector> dstVec;
  std::vector<MultiVector> scrapVec;
  for (int b = 0; b < blocks; b++) {
    dstVec.push_back(getBlock(b, dstCopy));
    scrapVec.push_back(getBlock(b, srcScrap_));
  }

  // run forward-substituion: run over each column
  //    From Heath pg. 65
  for (int b = 0; b < blocks; b++) {
    applyOp(invDiag_[b], scrapVec[b], dstVec[b]);

    // loop over each row
    for (int i = b + 1; i < blocks; i++) {
      LinearOp u_ib = getBlock(i, b, L_);
      if (u_ib != Teuchos::null) {
        applyOp(u_ib, dstVec[b], scrapVec[i], -1.0, 1.0);
      }
    }
  }

  // scale result by alpha
  if (beta != 0)
    update(alpha, dstCopy, beta, dst);  // dst = alpha * dstCopy + beta * dst
  else if (alpha != 1.0)
    scale(alpha, dst);  // dst = alpha * dst
}

void BlockLowerTriInverseOp::describe(Teuchos::FancyOStream& out_arg,
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
           << ",rows=" << blockRowCount(L_) << ",cols=" << blockColCount(L_) << "}\n";
      {
        OSTab tab2(out);
        *out << "[L Operator] = ";
        *out << Teuchos::describe(*L_, verbLevel);
      }
      {
        OSTab tab2(out);
        *out << "[invDiag Operators]:\n";
        tab.incrTab();
        for (int i = 0; i < blockRowCount(L_); i++) {
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
