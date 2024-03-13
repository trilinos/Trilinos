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
