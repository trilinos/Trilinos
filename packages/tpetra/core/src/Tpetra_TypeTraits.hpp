// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_TYPE_TRAITS_HPP
#define TPETRA_TYPE_TRAITS_HPP

/// \file Tpetra_TypeTraits.hpp

#include <type_traits>

namespace Tpetra {

/*! \brief is equivalent to std::true_type if memcpy can be called on this type.

    true if the type is trivially copyable
    true if the type defines an alias tpetra_memcpy_opt_in = std::true_type.
    This may be useful for types which are not trivially-copyable according to
   C++, but the implemented knows they are actually trivially copyable. Such a
   type S could add:

    struct S {
        ...
        using tpetra_memcpy_opt_in = std::true_type;
        ...
    }

    And tpetra will be willing to call memcpy rather than std::copy on that
   type.

    Alternatively, the implementer could provide an overload of okay_to_memcpy
   for the type in question

*/
template <typename T, typename Enable = void>
struct okay_to_memcpy : public std::false_type {};

template <typename T>
struct okay_to_memcpy<
    T, typename std::enable_if<std::is_trivially_copyable_v<T>>::type>
    : public std::true_type {};

// T::tpetra_memcpy_opt_in is defined as std::true_type.
template <typename T>
struct okay_to_memcpy<
    T, typename std::enable_if<std::is_same_v<typename T::tpetra_memcpy_opt_in,
                                              std::true_type>>::type>
    : public std::true_type {};

template <typename T>
constexpr bool okay_to_memcpy_v = okay_to_memcpy<T>::value;

} // namespace Tpetra

#endif /* TPETRA_TYPE_TRAITS */
