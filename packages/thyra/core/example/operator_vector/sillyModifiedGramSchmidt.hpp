// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
#define THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace Thyra {

/** \brief Silly little implementation of the modified Gram-Schmidt algorithm
 * to compute a QR factorization A=Q*R of a multi-vector A.
 *
 * \param V [in/out] On input, V contains the columns of A to compute the QR
 * factorization.  On output, V contains the columns of Q.
 *
 * \param R [out] On output, contains the upper triangular matrix R.
 *
 * \ingroup Thyra_Op_Vec_examples_cg_grp
 */
template<class Scalar>
void sillyModifiedGramSchmidt(
  const Ptr<MultiVectorBase<Scalar> > &V,
  const Ptr<RCP<MultiVectorBase<Scalar> > > &R_out
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST; using Teuchos::as;
  const int n = V->domain()->dim();
  *R_out = createMembers(V->domain(), V->domain());
  DetachedMultiVectorView<Scalar> R(*(*R_out));
  for (int k = 0; k < n; ++k) {
    R(k,k) = norm(*V->col(k));
    Vt_S(V->col(k).ptr(), ST::one()/R(k,k));
    for (int j = k+1; j < n; ++j) {
      R(k,j) = scalarProd(*V->col(k), *V->col(j));
      Vp_StV(V->col(j).ptr(), -R(k,j), *V->col(k));
    }
  }
} // end sillyModifiedGramSchmidt

} // namespace Thyra

#endif // THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
