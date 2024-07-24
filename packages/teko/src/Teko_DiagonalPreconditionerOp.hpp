// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_DiagonalPreconditionerOp_hpp__
#define __Teko_DiagonalPreconditionerOp_hpp__
#include "Teko_Utilities.hpp"

#include "Teko_ImplicitLinearOp.hpp"

class EpetraExt_PointToBlockDiagPermute;

namespace Teko {

class DiagonalPreconditionerOp : public ImplicitLinearOp {
 public:
  /** @brief Constuctor */
  DiagonalPreconditionerOp(Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP,
                           const VectorSpace range, const VectorSpace domain);

  /** @brief Range space of this operator */
  virtual VectorSpace range() const { return range_; }

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const { return domain_; }

  /** @brief Apply the preconditioner
   *
   * The <code>apply</code> function takes one vector as input
   * and applies a linear preconditioner. The result
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

  Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> get_BDP() const { return BDP_; }

 private:
  Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP_;
  const VectorSpace range_, domain_;
};

}  // end namespace Teko

#endif
