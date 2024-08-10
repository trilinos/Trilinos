// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_FUNCTIONBINDINGS_HPP
#define ROL_FUNCTIONBINDINGS_HPP

#include <functional>
#include "ROL_Vector.hpp"

namespace ROL {

namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
using f_update_t = function<void( const Vector<Real>& )>;

template<typename Real>
using f_scalar_t = function<Real( const Vector<Real>& )>;

template<typename Real>
using f_vector_t = function<void( Vector<Real>&, const Vector<Real>& )>;

template<typename Real>
using f_dderiv_t = function<void( Vector<Real>&, const Vector<Real>&, const Vector<Real>& )>;

template<typename Real>
using f_solve_t = function<void( Vector<Real> &, Vector<Real> & )>;

template<typename Real>
inline f_vector_t<Real> fix_direction( f_dderiv_t<Real>& f, const Vector<Real>& v ) {
  return bind( f, ph::_1, cref(v), ph::_2 );
}

template<typename Real>
inline f_vector_t<Real> fix_position( f_dderiv_t<Real>& f, const Vector<Real>& x ) {
  return bind( f, ph::_1, ph::_2, cref(x) );
}

} // namespace details 

template<typename Real> using f_update_t = details::f_update_t<Real>;
template<typename Real> using f_scalar_t = details::f_scalar_t<Real>;
template<typename Real> using f_vector_t = details::f_vector_t<Real>;
template<typename Real> using f_dderiv_t = details::f_dderiv_t<Real>;
template<typename Real> using f_solve_t  = details::f_solve_t<Real>;

template<typename Real>
inline f_vector_t<Real> fix_direction( f_dderiv_t<Real>& f, const Vector<Real>& v ) {
  return details::fix_direction(f,v);
}

template<typename Real>
inline f_vector_t<Real> fix_position( f_dderiv_t<Real>& f, const Vector<Real>& x ) {
  return details::fix_position(f,x);
}


} // namespace ROL


#endif // ROL_FUNCTIONBINDINGS_HPP

