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

#ifndef THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_DECL_HPP
#define THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_DECL_HPP

#include "Thyra_VectorSpaceConverterBase.hpp"

namespace Thyra {

/** \brief Node base class for converting serial multi-vectors (and vectors)
 * from one scalar type to another.
 *
 * This node base class defines the function <tt>convert()</tt> for all serial
 * vectors.  A concrete subclass is created by deriving from this interface
 * and then defining the function <tt>createVectorSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_serial_adapters_grp
 */
template<class ScalarFrom, class ScalarTo>
class SerialVectorSpaceConverterBase : virtual public VectorSpaceConverterBase<ScalarFrom,ScalarTo> {
public:

  /** @name Overridden from VectorSpaceConverterBase */
  //@{

  /** \brief . */
  void convert(
    const MultiVectorBase<ScalarFrom>    &mv_from
    ,MultiVectorBase<ScalarTo>           *mv_to
    ) const;

  //@}
  
};

} // namespace Thyra

#endif // THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_DECL_HPP
