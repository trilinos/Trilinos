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

#ifndef THYRA_VECTOR_SPACE_CONVERTED_BASE_HPP
#define THYRA_VECTOR_SPACE_CONVERTED_BASE_HPP

#include "Thyra_OperatorVectorTypes.hpp"


namespace Thyra {


/** \brief Base interface for a factory that converts vector space types and
 * vectors and multi-vectors from one scalar type to another.
 *
 * This abstract interface is templated on two scalar types,
 * <tt>ScalarFrom</tt> and <tt>ScalarTo</tt>, which allow a client to create
 * new vector spaces and convert vectors and multi-vectors from one scalar
 * type to another.  This is a critical capability for some algorithms with
 * mixed scalar types.
 *
 * One of these concrete subclasses should be provided for every concrete
 * <tt>VectorSpaceBase</tt> subclass.
 *
 * Note that this interface works just fine even if <tt>ScalarTo</tt> is the
 * same type as <tt>ScalarFrom</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class ScalarFrom, class ScalarTo>
class VectorSpaceConverterBase {
public:

  /** \brief . */
  virtual ~VectorSpaceConverterBase() {}

  /** \brief Create a new vector space with scalar type <tt>ScalarTo</tt>
   * given an existing vector space with scalar type <tt>ScalarFrom</tt>.
   */
  virtual Teuchos::RCP<const VectorSpaceBase<ScalarTo> >
  createVectorSpaceTo(
    const VectorSpaceBase<ScalarFrom>    &vecSpc
    ) const = 0;

  /** \brief Create a new vector space with scalar type <tt>ScalarFrom</tt>
   * given an existing vector space with scalar type <tt>ScalarTo</tt>.
   */
  virtual Teuchos::RCP<const VectorSpaceBase<ScalarFrom> >
  createVectorSpaceFrom(
    const VectorSpaceBase<ScalarTo>    &vecSpc
    ) const = 0;

  /** \brief Copy from a multi-vector (or vector) with one scalar type to
   * another multi-vector (or vector) with another scalar type.
   */
  virtual void convert(
    const MultiVectorBase<ScalarFrom>    &mv_from
    ,MultiVectorBase<ScalarTo>           *mv_to
    ) const = 0;
  
};


} // namespace Thyra


#endif // THYRA_VECTOR_SPACE_CONVERTED_BASE_HPP
