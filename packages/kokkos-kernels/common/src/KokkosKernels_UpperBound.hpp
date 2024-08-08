//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNELS_UPPERBOUND_HPP
#define _KOKKOSKERNELS_UPPERBOUND_HPP

/*! \file KokkosKernels_UpperBound.hpp
   Define thread and team-collaborative upper-bound search

   Upper-bound search takes a Kokkos::View, a search value, and a binary
   predicate.
   It returns an index to the first element of the view such that pred(value,
   element) is true

   This is implemented by calling lower_bound functions with inverted and
   reflected predicates, i.e. upper_bound(view, val, pred) = lower_bound(value,
   val, Inv(Refl(pred)));

   Examples:
   \verbatim
   value  = 3
   view   = {0,1,2,3,4}
          = {f,f,f,f,t}
   result =          4

   value  = -1
   view   = {0,1,2,3,4}
          = {t,t,t,t,t}
   result =  0

   value  = 5
   view   = {0,1,2,3,4}
          = {f,f,f,f,f}
   result =            5

   value  = 1
   view   = {0,1,1,1,2}
          = {f,f,f,f,t}
   result =          4
   \endverbatim

   Contrast with lower-bound, which returns first index for which pred(element,
   value) is false
 */

#include "KokkosKernels_LowerBound.hpp"

namespace KokkosKernels {

/*! \brief single-thread upper-bound search

    \tparam ViewLike A Kokkos::View or KokkosKernels::Impl::Iota
    \tparam Pred a binary predicate function
    \param view the view to search
    \param value the value to search for
    \param pred a binary predicate function
    \returns index of first element in view where pred(value,element) is true,
    or view.size if no such element exists
*/
template <typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type> >
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type upper_bound_thread(
    const ViewLike &view, const typename ViewLike::non_const_value_type &value, Pred pred = Pred()) {
  return lower_bound_thread(view, value, Neg(Refl(pred)));
}

/*! \brief team-collaborative upper-bound search

    \tparam ViewLike A Kokkos::View or KokkosKernels::Impl::Iota
    \tparam Pred a binary predicate function
    \param view the view to search
    \param value the value to search for
    \param pred a binary predicate function
    \returns index of first element in view where pred(value,element) is true,
    or view.size if no such element exists
*/
template <typename TeamMember, typename ViewLike, typename Pred = LT<typename ViewLike::non_const_value_type> >
KOKKOS_INLINE_FUNCTION typename ViewLike::size_type upper_bound_team(
    const TeamMember &handle, const ViewLike &view, const typename ViewLike::non_const_value_type &value,
    Pred pred = Pred()) {
  return lower_bound_team(handle, view, value, Neg(Refl(pred)));
}

}  // namespace KokkosKernels

#endif  // _KOKKOSKERNELS_UPPERBOUND_HPP