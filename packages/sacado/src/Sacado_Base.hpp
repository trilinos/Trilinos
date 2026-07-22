// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
