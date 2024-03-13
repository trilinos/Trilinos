/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#ifndef __Teko_BlockUpperTriInverseOp_hpp__
#define __Teko_BlockUpperTriInverseOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

/** \brief This linear operator computes the inverse
 *        of a upper triangular matrix.
 *
 * This linear operator computes the inverse
 * of an upper triangular matrix. This requires,
 * the upper triangular blocks, as well as the
 * inverse of the operators on the diagonal.
 */
class BlockUpperTriInverseOp : public BlockImplicitLinearOp {
 public:
  /** \brief This constructor explicitly takes an upper triangular matrix
   *        and inverse diagonal operators and builds a back substitution operator.
   *
   * This constructor explicitly takes an upper triangular matrix
   * and inverse diagonal operators and builds a back substitution operator.
   *
   * @param[in] U Upper triangular matrix object
   * @param[in] invDiag Vector containing the inverse of the diagonal blocks
   */
  BlockUpperTriInverseOp(BlockedLinearOp &U, const std::vector<LinearOp> &invDiag);

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

  /** @brief Perform a matrix vector multiply with this implicitly
   * defined blocked operator.
   *
   * The <code>apply</code> function takes one vector as input
   * and applies a linear operator. The result
   * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
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
  const BlockedLinearOp U_;        ///< operator \f$ U \f$
  std::vector<LinearOp> invDiag_;  ///< (Approximate) Inverses of the diagonal operators

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productDomain_;  ///< Domain vector space.

  // scratch space...so we don't have to reallocate
  mutable BlockedMultiVector srcScrap_;
  mutable BlockedMultiVector dstScrap_;

 private:
  // hide me!
  BlockUpperTriInverseOp();
  BlockUpperTriInverseOp(const BlockUpperTriInverseOp &);
};

inline LinearOp createBlockUpperTriInverseOp(BlockedLinearOp &U,
                                             const std::vector<LinearOp> &invDiag) {
  return Teuchos::rcp(new BlockUpperTriInverseOp(U, invDiag));
}

inline LinearOp createBlockUpperTriInverseOp(BlockedLinearOp &U,
                                             const std::vector<LinearOp> &invDiag,
                                             const std::string &str) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > result =
      Teuchos::rcp(new BlockUpperTriInverseOp(U, invDiag));
  result->setObjectLabel(str);

  return result;
}

}  // end namespace Teko

#endif
