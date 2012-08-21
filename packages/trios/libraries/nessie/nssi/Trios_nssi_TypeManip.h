/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// This code accompanies the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any
//     purpose is hereby granted without fee, provided that the above copyright
//     notice appear in all copies and that both that copyright notice and this
//     permission notice appear in supporting documentation.
// The author or Addison-Welsey Longman make no representations about the
//     suitability of this software for any purpose. It is provided "as is"
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

// Last update: November 22, 2001

#ifndef _TRIOS_NSSI_TYPEMANIP_H_
#define _TRIOS_NSSI_TYPEMANIP_H_

namespace Nessie
{
////////////////////////////////////////////////////////////////////////////////
// class template Int2Type
// Converts each integral constant into a unique type
// Invocation: Int2Type<v> where v is a compile-time constant integral
// Defines 'value', an enum that evaluates to v
////////////////////////////////////////////////////////////////////////////////

template <int v>
struct Int2Type
{
        enum { value = v };
};

////////////////////////////////////////////////////////////////////////////////
// class template Type2Type
// Converts each type into a unique, insipid type
// Invocation Type2Type<T> where T is a type
// Defines the type OriginalType which maps back to T
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Type2Type
{
        typedef T OriginalType;
};

////////////////////////////////////////////////////////////////////////////////
// class template Select
// Selects one of two types based upon a boolean constant
// Invocation: Select<flag, T, U>::Result
// where:
// flag is a compile-time boolean constant
// T and U are types
// Result evaluates to T if flag is true, and to U otherwise.
////////////////////////////////////////////////////////////////////////////////

template <bool flag, typename T, typename U>
struct Select
{
        typedef T Result;
};
template <typename T, typename U>
struct Select<false, T, U>
{
        typedef U Result;
};

////////////////////////////////////////////////////////////////////////////////
// class template IsSameType
// Return true iff two given types are the same
// Invocation: SameType<T, U>::value
// where:
// T and U are types
// Result evaluates to true iff U == T (types equal)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
struct IsSameType
{
        enum { value = false };
};

template <typename T>
struct IsSameType<T,T>
{
        enum { value = true };
};

////////////////////////////////////////////////////////////////////////////////
// Helper types Small and Big - guarantee that sizeof(Small) < sizeof(Big)
////////////////////////////////////////////////////////////////////////////////

namespace Private
{
template <class T, class U>
struct ConversionHelper
{
        typedef char Small;
        struct Big { char dummy[2]; };
        static Big   Test(...);
        static Small Test(U);
        static T MakeT();
};
}

////////////////////////////////////////////////////////////////////////////////
// class template Conversion
// Figures out the conversion relationships between two types
// Invocations (T and U are types):
// a) Conversion<T, U>::exists
// returns (at compile time) true if there is an implicit conversion from T
// to U (example: Derived to Base)
// b) Conversion<T, U>::exists2Way
// returns (at compile time) true if there are both conversions from T
// to U and from U to T (example: int to char and back)
// c) Conversion<T, U>::sameType
// returns (at compile time) true if T and U represent the same type
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

template <class T, class U>
struct Conversion
{
        typedef Private::ConversionHelper<T, U> H;
#ifndef __MWERKS__
        enum { exists = sizeof(typename H::Small) == sizeof((H::Test(H::MakeT()))) };
#else
        enum { exists = false };
#endif
        enum { exists2Way = exists && Conversion<U, T>::exists };
        enum { sameType = false };
};

template <class T>
struct Conversion<T, T>
{
        enum { exists = 1, exists2Way = 1, sameType = 1 };
};

template <class T>
struct Conversion<void, T>
{
        enum { exists = 0, exists2Way = 0, sameType = 0 };
};

template <class T>
struct Conversion<T, void>
{
        enum { exists = 0, exists2Way = 0, sameType = 0 };
};

template <>
struct Conversion<void, void>
{
    public:
        enum { exists = 1, exists2Way = 1, sameType = 1 };
};

////////////////////////////////////////////////////////////////////////////////
// class template SuperSubclass
// Invocation: SuperSubclass<B, D>::value where B and D are types.
// Returns true if B is a public base of D, or if B and D are aliases of the
// same type.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

template <class T, class U>
struct SuperSubclass
{
        enum { value = (::Nessie::Conversion<const volatile U*, const volatile T*>::exists &&
                !::Nessie::Conversion<const volatile T*, const volatile void*>::sameType) };
};

////////////////////////////////////////////////////////////////////////////////
// class template SuperSubclassStrict
// Invocation: SuperSubclassStrict<B, D>::value where B and D are types.
// Returns true if B is a public base of D.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

template<class T,class U>
struct SuperSubclassStrict
{
        enum { value = (::Nessie::Conversion<const volatile U*, const volatile T*>::exists &&
                !::Nessie::Conversion<const volatile T*, const volatile void*>::sameType &&
                !::Nessie::Conversion<const volatile T*, const volatile U*>::sameType) };
};

////////////////////////////////////////////////////////////////////////////////
// class template SizeOf
////////////////////////////////////////////////////////////////////////////////
template <class T>
class SizeOf
{
    public:
        static char check(Type2Type<void>);
        static char check(Type2Type<const void>);
        static char check(Type2Type<volatile void>);
        static char check(Type2Type<const volatile void>);
        static T check(...);

    public:
        enum {value = sizeof(check(Type2Type<T>()))};
};

}   // namespace Nessie

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS
// Invocation: SUPERSUBCLASS(B, D) where B and D are types.
// Returns true if B is a public base of D, or if B and D are aliases of the
// same type.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
// Deprecated: Use SuperSubclass class template instead.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS(T, U) \
        ::Nessie::SuperSubclass<T,U>::value

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS_STRICT
// Invocation: SUPERSUBCLASS(B, D) where B and D are types.
// Returns true if B is a public base of D.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
// Deprecated: Use SuperSubclassStrict class template instead.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS_STRICT(T, U) \
        ::Nessie::SuperSubclassStrict<T,U>::value

////////////////////////////////////////////////////////////////////////////////
// Change log:
// June 20, 2001: ported by Nick Thurn to gcc 2.95.3. Kudos, Nick!!!
// November 22, 2001: minor change to support porting to boost
// November 22, 2001: fixed bug in Conversion<void, T>
//      (credit due to Brad Town)
// November 23, 2001: (well it's 12:01 am) fixed bug in SUPERSUBCLASS - added
//      the volatile qualifier to be 100% politically correct
// September 16, 2002: Changed "const volatile" to "const volatile *", to enable
//     conversion to succeed. Done earlier by MKH.
//     Added SuperSubclass and SuperSubclassStrict templates. The corresponding
//     macros are deprecated.
//     Added extra parenthesis in sizeof in Conversion, to disambiguate function
//     call from function declaration. T.S.
////////////////////////////////////////////////////////////////////////////////

#endif // _TRIOS_NSSI_TYPEMANIP_H_
