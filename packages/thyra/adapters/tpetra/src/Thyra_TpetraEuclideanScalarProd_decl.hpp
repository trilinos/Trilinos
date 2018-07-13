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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
#define THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP


#include "Thyra_EuclideanScalarProd.hpp"
#include "Tpetra_MultiVector.hpp"


namespace Thyra {

/** \brief Extends concrete implementation of a Euclidean scalar product for
 * specifically Tpetra vectors/multivectors.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraEuclideanScalarProd : public EuclideanScalarProd<Scalar> {
protected:
  
  /** @name Overridden from EuclideanScalarProd */
  //@{

  /** \brief If X and Y are both Tpetra wrappers, computes the pair-wise
   * scalar products directly with Tpetra calls. Otherwise, this defers to
   * the base class implementaiton, which computes the result with an RTOp.
   */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X,
    const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar>& scalarProds
    ) const;

  //@}

private:

  /** /brief . */
  Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const;

};


/** \brief Nonmember constructor for TpetraEuclideanScalarProd.
 *
 * \relates TpetraEuclideanScalarProd
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<const TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraEuclideanScalarProd()
{
  return  Teuchos::rcp(new TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
}


} // end namespace Thyra


#endif  // THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
