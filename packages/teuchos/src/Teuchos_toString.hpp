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

#ifndef TEUCHOS_TO_STRING_HPP
#define TEUCHOS_TO_STRING_HPP

#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


/** \brief Default traits class for converting objects into strings.
 *
 * NOTE: This default implementation relies on opeator<<(std::ostream&, ...) 
 * being defined for the data type T.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
class ToStringTraits {
public:
  static std::string toString( const T &t )
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    }
};


/** \brief Utility function for returning a pretty string representation of
 * a object of type T.
 *
 * NOTE: This helper function simply returns ToStringTraits<T>::toString(t)
 * and the right way to speicalize the behavior is to specialize
 * ToStringTraits.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
inline
std::string toString(const T& t)
{
  return ToStringTraits<T>::toString(t);
}


/** \brief Specialization for bool. */
template<>
class ToStringTraits<bool> {
public:
  static std::string toString( const bool &t )
    {
      if (t)
        return "true";
      return "false";
    }
};


/** \brief Specialization for std::string. */
template<>
class ToStringTraits<std::string> {
public:
  static std::string toString( const std::string &t )
    {
      return t;
    }
};


} // end namespace Teuchos


#endif // TEUCHOS_TO_STRING_HPP
