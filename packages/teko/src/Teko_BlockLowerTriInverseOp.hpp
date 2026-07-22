// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockLowerTriInverseOp_hpp__
#define __Teko_BlockLowerTriInverseOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

/** \brief This linear operator computes the inverse
 *        of a lower triangular matrix.
 *
 * This linear operator computes the inverse
 * of a lower triangular matrix. This requires,
 * the lower triangular blocks, as well as the
 * inverse of the operators on the diagonal.
 */
class BlockLowerTriInverseOp : public BlockImplicitLinearOp {
 public:
  using BlockImplicitLinearOp::implicitApply;

  /** \brief This constructor explicitly takes a lower triangular matrix
   *        and inverse diagonal operators and builds a forward substitution operator.
   *
   * This constructor explicitly takes a lower triangular matrix
   * and inverse diagonal operators and builds a forward substitution operator.
   *
   * @param[in] L Lower triangular matrix object
   * @param[in] invDiag Vector containing the inverse of the diagonal blocks
   */
  BlockLowerTriInverseOp(BlockedLinearOp &L, const std::vector<LinearOp> &invDiag);

  //! @name Inherited methods from Thyra::LinearOpBase
  //@{

  /** @brief Range space of this operator */
  virtual VectorSpace range() const { return productRange_; }

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const { return productDomain_; }

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
  virtual void implicitApply(const BlockedMultiVector &x, BlockedMultiVector &y,
                             const double alpha = 1.0, const double beta = 0.0) const;
  //@}

  virtual void describe(Teuchos::FancyOStream &out_arg,
                        const Teuchos::EVerbosityLevel verbLevel) const;

 protected:
  // fundamental operators to use
  const BlockedLinearOp L_;        ///< operator \f$ L \f$
  std::vector<LinearOp> invDiag_;  ///< (Approximate) Inverses of the diagonal operators

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productDomain_;  ///< Domain vector space.

  // scratch space...so we don't have to reallocate
  mutable BlockedMultiVector srcScrap_;
  mutable BlockedMultiVector dstScrap_;
  mutable bool allocated = false;

 private:
  // hide me!
  BlockLowerTriInverseOp();
  BlockLowerTriInverseOp(const BlockLowerTriInverseOp &);
};

inline LinearOp createBlockLowerTriInverseOp(BlockedLinearOp &U,
                                             const std::vector<LinearOp> &invDiag) {
  return Teuchos::rcp(new BlockLowerTriInverseOp(U, invDiag));
}

inline LinearOp createBlockLowerTriInverseOp(BlockedLinearOp &L,
                                             const std::vector<LinearOp> &invDiag,
                                             const std::string &str) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > result =
      Teuchos::rcp(new BlockLowerTriInverseOp(L, invDiag));
  result->setObjectLabel(str);

  return result;
}

}  // end namespace Teko

#endif
