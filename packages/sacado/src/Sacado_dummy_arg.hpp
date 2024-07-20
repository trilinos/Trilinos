// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
