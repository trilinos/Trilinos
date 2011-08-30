// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
 * \file   Amesos2_Meta.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Mon Jul 18 15:28:52 2011
 *
 * \brief  Provides some simple meta-programming utilities for Amesos2.
 *
 * The meta-programming utilities provided here might not be nearly as
 * flexible and complete as those provided in the Boost
 * Meta-Programming Library, but they are sufficient for the needs of
 * Amesos2.
 */


#ifndef AMESOS2_META_HPP
#define AMESOS2_META_HPP

#include "Amesos2_config.h"


namespace Amesos2 {

  namespace Meta {

    /**
     * \internal
     * \defgroup amesos2_meta Amesos2 Meta-programming Functions
     * @{
     */

    template <class T, T val>
    struct integral_constant
    {
      typedef integral_constant<T, val>  type;
      typedef T                          value_type;
      static const T value;
    };

    /* Some compilers support initializing static const members alongside the
     * definition, but others do not, so we go we the safe method of external
     * initialization.
     */
    template <class T, T val>
    const T integral_constant<T,val>::value = val;

    typedef integral_constant<bool, true>  true_type;
    typedef integral_constant<bool, false> false_type;

    ////////////////////////////////////////
    // Testing the same'ness of two types //
    ////////////////////////////////////////

    /**
     * \brief test same-ness of two types
     */
    template <typename, typename>
    struct is_same : public false_type
    {};

    template <typename T>
    struct is_same<T,T> : public true_type
    {};



    //////////////////////////////////////////
    // Meta-functions for boolean operators //
    //////////////////////////////////////////

    /* Must define these with the '_' suffix because they are
     * otherwise keywords in C++
     */

    template <bool b1, bool b2>
    struct or_ : public false_type {};

    template <bool b>
    struct or_<true,b> : public true_type {};

    template <bool b>
    struct or_<b,true> : public true_type {};


    template <bool b1, bool b2>
    struct and_ : public false_type {};

    template <>
    struct and_<true,true> : public true_type {};


    template <bool b>
    struct not_ {};

    template <>
    struct not_<true> : false_type {};

    template <>
    struct not_<false> : true_type {};


    //////////////////////////////////////
    // Evaluating to a conditional type //
    //////////////////////////////////////

    template <bool B, typename T1, typename T2>
    struct if_then_else {};

    template <typename T1, typename T2>
    struct if_then_else<true, T1, T2> {
      typedef T1 type;
    };

    template <typename T1, typename T2>
    struct if_then_else<false, T1, T2> {
      typedef T2 type;
    };


    ////////////////////////////////////////////
    // A meta-programming type-list structure //
    ////////////////////////////////////////////

    struct nil_t {};                // to denote an empty list

    template <typename Head, typename Tail>
    struct type_list {
      typedef type_list<Head,Tail> type;
      typedef Head head;
      typedef Tail tail;
    };

    /**
     * \brief Utility meta-function for creating a type-list
     *
     * Example:
     *
     * \code
     * make_list3<float,double,quad>
     * \encode
     *
     * is a list of three types: \c float , \c double , and \c quad .
     */
    template <typename T1>
    struct make_list1
      : type_list<T1,nil_t>
    { };

    template <typename T1, typename T2>
    struct make_list2
      : type_list<T1,type_list<T2,nil_t> >
    { };

    template <typename T1, typename T2, typename T3>
    struct make_list3
      : type_list<T1, type_list<T2, type_list<T3,nil_t> > >
    { };

    template <typename T1, typename T2, typename T3, typename T4>
    struct make_list4
      : type_list<T1, type_list<T2, type_list<T3, type_list<T4,nil_t> > > >
    { };

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    struct make_list5
      : type_list<T1, type_list<T2, type_list<T3, type_list<T4, type_list<T5,nil_t> > > > >
    { };

    template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    struct make_list6
      : type_list<T1, type_list<T2, type_list<T3, type_list<T4, type_list<T5, type_list<T6,nil_t> > > > > >
    { };

    /* More declarations for larger type lists may be added if necessary */


    /**
     * \brief A utility meta-function to determine whether a given
     * type is found within a type-list
     *
     * \code
     * typedef make_list4<float,double,int,long> t_list;
     * if( type_list_contains<t_list,bool>::value ){
     *   // dead branch
     * } else {
     *   // This will always execute
     * }
     */
    template <typename list, typename elem>
    struct type_list_contains {
      static const bool value;
    };

    // Base recursive case
    template <typename elem>
    struct type_list_contains<nil_t,elem> {
      static const bool value = false;
    };

    template <typename list, typename elem>
    const bool type_list_contains<list,elem>::value
    = if_then_else<is_same<typename list::head, elem>::value,
                   true_type,
                   type_list_contains<typename list::tail,elem> >::type::value;

    /** @} */

  } // end namespace Meta

} // end namespace Amesos2

#endif  // AMESOS2_META_HPP
