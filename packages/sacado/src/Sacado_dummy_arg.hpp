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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_DUMMY_ARG_HPP
#define SACADO_DUMMY_ARG_HPP

namespace Sacado {

  //! A dummy argument that can be converted to any scalar type
  template <class T> struct dummy_arg { 
    operator T() const { return T(0.0); } 
  };

  //! A meta-function that defines U as its type
  template <class T, class U> struct dummy { 
    typedef U type; 
  };

  //! Specialization to provide a dummy argument when types are the same
  template <class T> struct dummy<T,T> { 
    typedef dummy_arg<T> type; 
  };

} // namespace Sacado

#endif // SACADO_DUMMY_ARG_HPP
