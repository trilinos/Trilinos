// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef L1_UTILITIES_HPP
#define L1_UTILITIES_HPP

#include <tuple>
//#include "ROL_BoundConstraint.hpp"
//#include "ROL_Constraint.hpp"
//#include "ROL_Objective.hpp"
#include "ROL_PartitionedVector.hpp"
//#include "ROL_Ptr.hpp"

namespace ROL {

//template<typename> class L1PenaltyObjective;
//template<typename> class L1PenaltyConstraint;

namespace detail { 

template<typename Real>
inline auto pv_cast( const Vector<Real>& x ) {
  return static_cast<const PartitionedVector<Real>&>(x);
}

template<typename Real>
inline auto pv_cast( Vector<Real>& x ) {
  return static_cast<PartitionedVector<Real>&>(x);
}

//template<template<typename> class Derived, typename Real>
//struct Base {};
//
//template<typename Real>
//struct Base<PartitionedVector,Real> { using type = Vector<Real>; };
//
//template<typename Real>
//struct Base<BoundConstraint_Partitioned,Real> { using type = BoundConstraint<Real>; };
//
//template<typename Real>
//struct Base<L1PenaltyConstraint,Real> { using type = Constraint<Real>; };
//
//template<typename Real>
//struct Base<L1PenaltyObjective,Real> { using type = Objective<Real>; };


} // namespace detail 

auto&& pv_cast = []( auto&&...args ) { return std::make_tuple(detail::pv_cast(args)...); };

//template<typename Real>
//inline auto up_cast( const Ptr<const PartitionedVector<Real>>& x ) { 
//  return staticPtrCast<const Vector<Real>>(x);
//}
//
//template<typename Real>
//inline auto up_cast( const Ptr<PartitionedVector<Real>>& x ) { 
//  return staticPtrCast<Vector<Real>>(x);
//}
//
//template<template<typename> class C, typename T>
//using base_t = typename detail::Base<C,T>::type;
//
//template<template<typename> class Derived, typename Real>
//inline auto up_cast( const Ptr<const Derived<Real>>& x ) { 
//  return staticPtrCast<const base_t<Derived,Real>>(x);
//}
//
//template<template<typename> class Derived, typename Real>
//inline auto up_cast( const Ptr<Derived<Real>>& x ) { 
//  return staticPtrCast<base_t<Derived,Real>>(x);
//}


//template<typename Real>
//inline auto up_cast( const Ptr<const PartitionedVector<Real>>& x ) { 
//  return staticPtrCast<const Vector<Real>>(x);
//}
//
//template<typename Real>
//inline auto up_cast( const Ptr<PartitionedVector<Real>>& x ) { 
//  return staticPtrCast<Vector<Real>>(x);
//}



} // namespace ROL



#endif //L1_UTILITIES_HPP

