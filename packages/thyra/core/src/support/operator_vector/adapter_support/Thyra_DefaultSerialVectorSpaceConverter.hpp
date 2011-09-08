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

#ifndef THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP
#define THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP

#include "Thyra_SerialVectorSpaceConverterBase.hpp"


namespace Thyra {


/** \brief Concrete subclass for a converter subclass for converting serial
 * multi-vectors and vectors.
 *
 * While this concrete subclass creates concrete vector spaces of type
 * <tt>DefaultSerialVectorSpace</tt>, it should be usable with any serial vector
 * space type and therefore this subclass is more general then it may appear
 * at first.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class ScalarFrom, class ScalarTo>
class DefaultSerialVectorSpaceConverter : virtual public SerialVectorSpaceConverterBase<ScalarFrom,ScalarTo> {
public:

  /** @name Overridden from VectorSpaceConverterBase */
  //@{

  /** \brief . */
  virtual
  Teuchos::RCP<const VectorSpaceBase<ScalarTo> >
  createVectorSpaceTo(
    const VectorSpaceBase<ScalarFrom> &vecSpc
    ) const;

  /** \brief . */
  virtual
  Teuchos::RCP<const VectorSpaceBase<ScalarFrom> >
  createVectorSpaceFrom(
    const VectorSpaceBase<ScalarTo> &vecSpc
    ) const;

  //@}
  
};


// Implementation


template<class ScalarFrom, class ScalarTo>
Teuchos::RCP<const VectorSpaceBase<ScalarTo> >
DefaultSerialVectorSpaceConverter<ScalarFrom,ScalarTo>::createVectorSpaceTo(
  const VectorSpaceBase<ScalarFrom>  &vecSpc
  ) const
{
  return defaultSpmdVectorSpace<ScalarTo>(vecSpc.dim());
}


template<class ScalarFrom, class ScalarTo>
Teuchos::RCP<const VectorSpaceBase<ScalarFrom> >
DefaultSerialVectorSpaceConverter<ScalarFrom,ScalarTo>::createVectorSpaceFrom(
  const VectorSpaceBase<ScalarTo>  &vecSpc
  ) const
{
  return defaultSpmdVectorSpace<ScalarFrom>(vecSpc.dim());
}


} // namespace Thyra


#endif // THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP
