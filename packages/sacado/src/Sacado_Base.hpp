// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_BASE_HPP
#define SACADO_BASE_HPP

namespace Sacado {

  //! Base class for Sacado types to control overload resolution
  /*!
   * Non-expression template-based Sacado AD types should be derived from this
   * class using the AD type as the template parameter (through CRTP) and
   * implement their overloads using this type instead of the AD type.  The
   * purpose of this is to control the visible overload set for the math
   * functions so Sacado's overloads only match when at least one of the
   * arguments is in fact a Sacado type (since most AD types enable implicit
   * conversions from their inner scalar type).
   */
  template <typename T>
  struct Base {
    typedef T derived_type;
    const derived_type& derived() const {
      return static_cast<const derived_type&>(*this);
    }
  };

}

#endif // SACADO_BASE_HPP
