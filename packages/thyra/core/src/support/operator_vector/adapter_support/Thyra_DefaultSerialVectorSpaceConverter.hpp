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
