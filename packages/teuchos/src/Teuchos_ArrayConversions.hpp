// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_ARRAY_CONVERSIONS_H
#define TEUCHOS_ARRAY_CONVERSIONS_H

/*! \file Teuchos_ArrayConversions.hpp
  \brief Templated conversions between Array<RCP<T> > and ArrayView<Ptr<T> >
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"


namespace Teuchos {


/** \brief Utility function to convert Array<RCP<T> > to Array<Ptr<T> >,
 * and Array<RCP<T> > to Array<Ptr<const T> >.
 *
 *  Templated on an array of pointer objects on input
 */

template<class T_out, class ArrayPtrT_in>
Array<Ptr<T_out> > arrayPtrConv(const ArrayPtrT_in &a_in)
{
  using Teuchos::as;
  Array<Ptr<T_out> > a_out(a_in.size());
  for (Teuchos_Ordinal i = 0; i < as<Teuchos_Ordinal>(a_in.size()); ++i) {
    a_out[i] = a_in[i].ptr();
  }
  return a_out;
}

} // namespace Teuchos

#endif // TEUCHOS_ARRAY_CONVERSIONS_H

