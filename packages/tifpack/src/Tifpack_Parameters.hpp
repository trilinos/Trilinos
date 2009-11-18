/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef TIFPACK_PARAMETERS_HPP
#define TIFPACK_PARAMETERS_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tifpack {

//! Fills a list which contains all the parameters possibly used by Tifpack.
void GetValidParameters(Teuchos::ParameterList& params);

//! Set a value from a ParameterList if a parameter with the specified name exists.
/** If the specified name does not name a parameter in the list, then 'value' is
  not referenced.
*/
template<typename T>
void GetParameter(const Teuchos::ParameterList& params, const std::string& name, T& value)
{
  if (params.isParameter(name)) {
    if (params.isType<T>(name)) {
      value = params.get<T>(name);
    }
  }
}

}//namespace Tifpack

#endif
