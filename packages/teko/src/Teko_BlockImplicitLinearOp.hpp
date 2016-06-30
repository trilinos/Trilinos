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
   virtual void implicitApply(const BlockedMultiVector & x, BlockedMultiVector & y,
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
   virtual void implicitApply(const Thyra::EOpTransp M_trans,
              const BlockedMultiVector & x, BlockedMultiVector & y,
              const double alpha = 1.0, const double beta = 0.0) const;

protected:

   //! Functions required by Thyra::LinearOpBase 
   //@{ 

  virtual bool opSupportedImpl(const Thyra::EOpTransp M_trans) const;

  virtual void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<double> & x,
    const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & y,
    const double alpha,
    const double beta
    ) const;
  
   //@}
};

} // end namespace Teko

#endif
