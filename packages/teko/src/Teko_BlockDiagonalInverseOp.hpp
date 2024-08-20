// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockDiagonalInverseOp_hpp__
#define __Teko_BlockDiagonalInverseOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

/** \brief This linear operator computes the inverse
 *        of a diagonal matrix.
 *
 * This linear operator computes the inverse
 * of a diagonal matrix. This requires
 * the inverse of the operators on the diagonal.
 */
class BlockDiagonalInverseOp : public BlockImplicitLinearOp {
 public:
  /** \brief This constructor explicitly takes a diagonal matrix
   *        and inverse diagonal operators and builds a back substitution operator.
   *
   * This constructor explicitly takes the
   * inverse diagonal operators and builds a back substitution operator.
   *
   * @param[in] A Matrix object
   * @param[in] invDiag Vector containing the inverse of the diagonal blocks
   */
  BlockDiagonalInverseOp(BlockedLinearOp &A, const std::vector<LinearOp> &invDiag);

  BlockDiagonalInverseOp()                               = delete;
  BlockDiagonalInverseOp(const BlockDiagonalInverseOp &) = delete;

  //! @name Inherited methods from Thyra::LinearOpBase
  //@{

  /** @brief Range space of this operator */
  virtual VectorSpace range() const { return productRange_; }

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const { return productDomain_; }

  /** @brief Perform a matrix vector multiply with this operator.
   *
   * The <code>apply</code> function takes one vector as input
   * and applies \f$ D^{-1} \f$. The result
   * is returned in \f$y\f$. If this operator is represented as \f$M\f$ then
   * \f$ y = \alpha M x + \beta y \f$.
   *
   * @param[in]     x
   * @param[in,out] y
   * @param[in]     alpha (default=1)
   * @param[in]     beta  (default=0)
   */
  virtual void implicitApply(const BlockedMultiVector &x, BlockedMultiVector &y,
                             const double alpha = 1.0, const double beta = 0.0) const;

  /** @brief Perform a matrix vector multiply with this implicitly
   * defined blocked operator.
   *
   * The <code>apply</code> function takes one vector as input
   * and applies a linear operator. The result
   * is returned in \f$y\f$. If this operator is represented as \f$M\f$ then
   * \f$ y = \alpha M x + \beta y \f$
   *
   * @param[in]     x
   * @param[in,out] y
   * @param[in]     alpha (default=1)
   * @param[in]     beta  (default=0)
   */
  virtual void implicitApply(const Thyra::EOpTransp M_trans, const BlockedMultiVector &x,
                             BlockedMultiVector &y, const double alpha = 1.0,
                             const double beta = 0.0) const;
  //@}

  virtual void describe(Teuchos::FancyOStream &out_arg,
                        const Teuchos::EVerbosityLevel verbLevel) const;

 protected:
  // fundamental operators to use
  std::vector<LinearOp> invDiag_;  ///< (Approximate) Inverses of the diagonal operators

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productDomain_;  ///< Domain vector space.

  // scratch space...so we don't have to reallocate
  mutable BlockedMultiVector srcScrap_;
  mutable BlockedMultiVector dstScrap_;
  mutable bool allocated = false;
};

inline LinearOp createBlockDiagonalInverseOp(BlockedLinearOp &A,
                                             const std::vector<LinearOp> &invDiag) {
  return Teuchos::rcp(new BlockDiagonalInverseOp(A, invDiag));
}

inline LinearOp createBlockDiagonalInverseOp(BlockedLinearOp &A,
                                             const std::vector<LinearOp> &invDiag,
                                             const std::string &str) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > result =
      Teuchos::rcp(new BlockDiagonalInverseOp(A, invDiag));
  result->setObjectLabel(str);

  return result;
}

}  // end namespace Teko

#endif
