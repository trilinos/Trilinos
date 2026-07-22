// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockImplicitLinearOp_hpp__
#define __Teko_BlockImplicitLinearOp_hpp__

#include "Teko_Utilities.hpp"

namespace Teko {

/** \brief A virtual class that simplifies the construction
 *        of custom operators.
 *
 * A virtual class that simplifies the construction
 * of custom operators. Good examples can be found in <code>LU2x2InverseOp</code>
 * and in <code>BlockUpperTriInverseOp</code>.
 */
class BlockImplicitLinearOp : public Thyra::LinearOpBase<double> {
 public:
  /** @brief Range space of this operator */
  virtual VectorSpace range() const = 0;

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const = 0;

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
  virtual void implicitApply(const BlockedMultiVector& x, BlockedMultiVector& y,
                             const double alpha = 1.0, const double beta = 0.0) const = 0;

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
  virtual void implicitApply(const Thyra::EOpTransp M_trans, const BlockedMultiVector& x,
                             BlockedMultiVector& y, const double alpha = 1.0,
                             const double beta = 0.0) const;

 protected:
  //! Functions required by Thyra::LinearOpBase
  //@{

  virtual bool opSupportedImpl(const Thyra::EOpTransp M_trans) const;

  virtual void applyImpl(const Thyra::EOpTransp M_trans, const Thyra::MultiVectorBase<double>& x,
                         const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& y, const double alpha,
                         const double beta) const;

  //@}
};

}  // end namespace Teko

#endif
