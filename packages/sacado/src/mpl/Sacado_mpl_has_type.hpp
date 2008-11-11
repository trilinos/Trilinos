// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_MPL_HAS_TYPE_HPP
#define SACADO_MPL_HAS_TYPE_HPP

namespace Sacado {

  namespace mpl {

    // Uses SFINAE to determine if a type has a nested typedef called "type",
    // i.e., whether the type is a metafunction

    typedef char NotFound;           // sizeof(NotFound) == 1
    struct Found { char x[2]; };     // sizeof(Found) == 2

    // Overload matches if and only if T::type exists
    template <class T> Found testHasType(typename T::type*);

    // Compiler prefers anything at all over ...
    template <class T> NotFound testHasType(...);

    template <class T> struct has_type {
      static const bool value = (sizeof(testHasType<T>(0)) == sizeof(Found));
    };

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_HAS_TYPE_HPP
