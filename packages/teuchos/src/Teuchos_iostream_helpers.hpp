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

#ifndef TEUCHOS_IOSTREAM_HELPERS_HPP
#define TEUCHOS_IOSTREAM_HELPERS_HPP


/*! \file Teuchos_VerbosityLevel.hpp
    \brief .
*/


#include "Teuchos_TestForException.hpp"


namespace Teuchos {


template<class EnumType>
std::istream& enumIstreamExtractor(std::istream& std_is, EnumType& enum_value)
{
  int intval;
  std_is >> intval;
  enum_value = static_cast<EnumType>(intval);
  return std_is;
}


} // namespace Teuchos


#define TEUCHOS_ENUM_INPUT_STERAM_OPERATOR(ENUMTYPE) \
inline \
std::istream& operator>>(std::istream& std_is, ENUMTYPE& enum_value) \
{ return Teuchos::enumIstreamExtractor(std_is, enum_value); }


#endif // TEUCHOS_IOSTREAM_HELPERS_HPP
