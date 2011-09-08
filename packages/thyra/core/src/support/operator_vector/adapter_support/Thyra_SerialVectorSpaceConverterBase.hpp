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

#ifndef THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP
#define THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP

#include "Thyra_VectorSpaceConverterBase.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"


namespace Thyra {


/** \brief Node base class for converting serial multi-vectors (and vectors)
 * from one scalar type to another.
 *
 * This node base class defines the function <tt>convert()</tt> for all serial
 * vectors.  A concrete subclass is created by deriving from this interface
 * and then defining the function <tt>createVectorSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template<class ScalarFrom, class ScalarTo>
class SerialVectorSpaceConverterBase
  : virtual public VectorSpaceConverterBase<ScalarFrom,ScalarTo>
{
public:

  /** @name Overridden from VectorSpaceConverterBase */
  //@{

  /** \brief . */
  virtual void convert(
    const MultiVectorBase<ScalarFrom> &mv_from,
    MultiVectorBase<ScalarTo> *mv_to
    ) const;

  //@}
  
};


// Implementation


template<class ScalarFrom, class ScalarTo>
void SerialVectorSpaceConverterBase<ScalarFrom,ScalarTo>::convert(
  const MultiVectorBase<ScalarFrom> &mv_from,
  MultiVectorBase<ScalarTo> *mv_to
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(mv_to==NULL);
#endif
  ConstDetachedMultiVectorView<ScalarFrom> emv_from(mv_from);
  DetachedMultiVectorView<ScalarTo> emv_to(*mv_to);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(emv_from.subDim() != emv_to.subDim());
  TEST_FOR_EXCEPT(emv_from.numSubCols() != emv_to.numSubCols());
#endif
  for( Ordinal j = 0; j < emv_from.numSubCols(); ++j ) {
    for( Ordinal i = 0; i < emv_from.subDim(); ++i ) {
      emv_to(i,j) = emv_from(i,j); // ToDo: Make this faster using optimized copy functions?
    }
  }
}


} // namespace Thyra


#endif // THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP
