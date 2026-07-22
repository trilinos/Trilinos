// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
