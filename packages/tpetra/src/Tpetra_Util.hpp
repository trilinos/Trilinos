// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_UTIL_HPP
#define TPETRA_UTIL_HPP

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 
#include <Teuchos_Utils.hpp>

namespace Tpetra {

  // efficientAddOrUpdate is taken from Scott Meyers' "Effective STL", Item 24.
  // if m already contains an entry with key k, use operator [].
  // if it doesn't, insert is used.
	template<typename MapType, typename KeyArgType, typename ValueArgType>
  typename MapType::iterator efficientAddOrUpdate(MapType& m, KeyArgType const& k, ValueArgType const& v) {
    typename MapType::iterator lb = m.lower_bound(k);
    if(lb != m.end() && !(m.key_comp()(k, lb->first))) {
      lb->second = v;
      return(lb);
    }
    else {
      typedef typename MapType::value_type MVT;
      return(m.insert(lb, MVT(k, v)));
    }
  }

  // type conversion functions
  template<typename OrdinalType, typename ScalarType>
  void ordinalToScalar(OrdinalType const& source, ScalarType& dest) {
    dest = static_cast<ScalarType>(source);
  }
  template<typename OrdinalType, typename ScalarType>
  void scalarToOrdinal(ScalarType const& source, OrdinalType& dest) {
    dest = static_cast<OrdinalType>(source);
  }

  // functions for converting types to strings
  // mainly used for doing output
  template <typename T>
  std::string toString(T const& x)
  {
    return(Teuchos::toString(x));
  }

template<typename T>
inline ostream& operator<<(ostream& os, std::vector<T> const& vector)
{
  os << "{";
  if(!vector.empty()) {
    typename std::vector<T>::const_iterator i = vector.begin();
    os << *i;
    i++;
    for(; i != vector.end(); i++)
      os << "," << *i;
  }
  os << "}";
  return(os);
}

} // namespace Tpetra


// this function works much the way Teuchos::Array::toString works.
// it allows std::vector to be used with an ostream.
// The contents of the vector are printed in the following format:
// "{4, 7, 18, 23, 6, 2}"
template<typename T>
inline ostream& operator<<(ostream& os, std::vector<T> const& vector)
{
  os << "{";
  if(!vector.empty()) {
    typename std::vector<T>::const_iterator i = vector.begin();
    os << *i;
    i++;
    for(; i != vector.end(); i++)
      os << "," << *i;
  }
  os << "}";
  return(os);
}

#endif // TPETRA_UTIL_HPP
