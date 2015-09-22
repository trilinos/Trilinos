/*
Copyright (C) 2005,2009-2010 Electronic Arts, Inc.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1.  Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
3.  Neither the name of Electronic Arts, Inc. ("EA") nor the names of
    its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY ELECTRONIC ARTS AND ITS CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ELECTRONIC ARTS OR ITS CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

///////////////////////////////////////////////////////////////////////////////
// EASTL/internal/type_fundamental.h
//
// Written and maintained by Paul Pedriana - 2005.
///////////////////////////////////////////////////////////////////////////////


#ifndef EASTL_INTERNAL_TYPE_FUNDAMENTAL_H
#define EASTL_INTERNAL_TYPE_FUNDAMENTAL_H


namespace eastl
{

    // The following properties or relations are defined here. If the given
    // item is missing then it simply hasn't been implemented, at least not yet.


    ///////////////////////////////////////////////////////////////////////
    // is_void
    //
    // is_void<T>::value == true if and only if T is one of the following types:
    //    [const][volatile] void
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_void : public false_type{};

    template <> struct is_void<void> : public true_type{};
    template <> struct is_void<void const> : public true_type{};
    template <> struct is_void<void volatile> : public true_type{};
    template <> struct is_void<void const volatile> : public true_type{};


    ///////////////////////////////////////////////////////////////////////
    // is_integral
    //
    // is_integral<T>::value == true if and only if T  is one of the following types:
    //    [const] [volatile] bool
    //    [const] [volatile] char
    //    [const] [volatile] signed char
    //    [const] [volatile] unsigned char
    //    [const] [volatile] wchar_t
    //    [const] [volatile] short
    //    [const] [volatile] int
    //    [const] [volatile] long
    //    [const] [volatile] long long
    //    [const] [volatile] unsigned short
    //    [const] [volatile] unsigned int
    //    [const] [volatile] unsigned long
    //    [const] [volatile] unsigned long long
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_integral : public false_type{};

    // To do: Need to define volatile and const volatile versions of these.
    template <> struct is_integral<unsigned char>            : public true_type{};
    template <> struct is_integral<const unsigned char>      : public true_type{};
    template <> struct is_integral<unsigned short>           : public true_type{};
    template <> struct is_integral<const unsigned short>     : public true_type{};
    template <> struct is_integral<unsigned int>             : public true_type{};
    template <> struct is_integral<const unsigned int>       : public true_type{};
    template <> struct is_integral<unsigned long>            : public true_type{};
    template <> struct is_integral<const unsigned long>      : public true_type{};
    template <> struct is_integral<unsigned long long>       : public true_type{};
    template <> struct is_integral<const unsigned long long> : public true_type{};

    template <> struct is_integral<signed char>              : public true_type{};
    template <> struct is_integral<const signed char>        : public true_type{};
    template <> struct is_integral<signed short>             : public true_type{};
    template <> struct is_integral<const signed short>       : public true_type{};
    template <> struct is_integral<signed int>               : public true_type{};
    template <> struct is_integral<const signed int>         : public true_type{};
    template <> struct is_integral<signed long>              : public true_type{};
    template <> struct is_integral<const signed long>        : public true_type{};
    template <> struct is_integral<signed long long>         : public true_type{};
    template <> struct is_integral<const signed long long>   : public true_type{};

    template <> struct is_integral<bool>            : public true_type{};
    template <> struct is_integral<const bool>      : public true_type{};
    template <> struct is_integral<char>            : public true_type{};
    template <> struct is_integral<const char>      : public true_type{};
    #ifndef EA_WCHAR_T_NON_NATIVE // If wchar_t is a native type instead of simply a define to an existing type...
        template <> struct is_integral<wchar_t>         : public true_type{};
        template <> struct is_integral<const wchar_t>   : public true_type{};
    #endif

    ///////////////////////////////////////////////////////////////////////
    // is_floating_point
    //
    // is_floating_point<T>::value == true if and only if T is one of the following types:
    //    [const] [volatile] float
    //    [const] [volatile] double
    //    [const] [volatile] long double
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_floating_point : public false_type{};

    // To do: Need to define volatile and const volatile versions of these.
    template <> struct is_floating_point<float>             : public true_type{};
    template <> struct is_floating_point<const float>       : public true_type{};
    template <> struct is_floating_point<double>            : public true_type{};
    template <> struct is_floating_point<const double>      : public true_type{};
    template <> struct is_floating_point<long double>       : public true_type{};
    template <> struct is_floating_point<const long double> : public true_type{};



    ///////////////////////////////////////////////////////////////////////
    // is_arithmetic
    //
    // is_arithmetic<T>::value == true if and only if:
    //    is_floating_point<T>::value == true, or
    //    is_integral<T>::value == true
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_arithmetic : public integral_constant<bool,
        is_integral<T>::value || is_floating_point<T>::value
    >{};


    ///////////////////////////////////////////////////////////////////////
    // is_fundamental
    //
    // is_fundamental<T>::value == true if and only if:
    //    is_floating_point<T>::value == true, or
    //    is_integral<T>::value == true, or
    //    is_void<T>::value == true
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_fundamental : public integral_constant<bool,
        is_void<T>::value || is_integral<T>::value || is_floating_point<T>::value
    >{};

} // namespace eastl


#endif // Header include guard
