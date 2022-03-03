// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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

