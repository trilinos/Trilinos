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

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {

	// this will contain utility functions, such as efficientAddOrUpdate

  // efficientAddOrUpdate is taken from Scott Meyers' "Effective STL", Item 24.
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

} // namespace Tpetra

#endif // TPETRA_UTIL_HPP
