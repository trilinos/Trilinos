// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
