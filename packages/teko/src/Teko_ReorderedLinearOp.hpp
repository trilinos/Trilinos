// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_ReorderedLinearOp_hpp__
#define __Teko_ReorderedLinearOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_ImplicitLinearOp.hpp"
#include "Teko_BlockedReordering.hpp"

namespace Teko {

/** \brief This class takes a blocked linear op and represents it
 *        in a flattened form.
 */
class ReorderedLinearOp : public ImplicitLinearOp {
 public:
  ReorderedLinearOp(const Teuchos::RCP<const BlockReorderManager> &mgr,
                    const Teuchos::RCP<Thyra::LinearOpBase<double> > &blockedOp);

  //! Get accessor to make reuse easier
  Teuchos::RCP<const BlockReorderManager> getReorderManager() const { return mgr_; }

  //! Get accessor to make reuse easier
  Teko::ModifiableLinearOp getBlockedOp() const { return blockedOp_; }

  /** @brief Range space of this operator */
  virtual VectorSpace range() const;

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const;

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
  virtual void implicitApply(const MultiVector &x, MultiVector &y, const double alpha = 1.0,
                             const double beta = 0.0) const;

  virtual void describe(Teuchos::FancyOStream &out_arg,
                        const Teuchos::EVerbosityLevel verbLevel) const;

 private:
  VectorSpace range_;
  VectorSpace domain_;
  Teuchos::RCP<const BlockReorderManager> mgr_;
  Teuchos::RCP<Thyra::LinearOpBase<double> > blockedOp_;
};

}  // end namespace Teko

#endif
