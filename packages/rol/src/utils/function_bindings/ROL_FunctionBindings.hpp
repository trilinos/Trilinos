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
#ifndef ROL_FUNCTIONBINDING_HPP
#define ROL_FUNCTIONBINDING_HPP

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
f_vector_t<Real> fix_direction( f_dderiv_t<Real>& f, const Vector<Real>& v ) {
  return bind( f, ph::_1, cref(v), ph::_2 );
}

template<typename Real>
f_vector_t<Real> fix_position( f_dderiv_t<Real>& f, const Vector<Real>& x ) {
  return bind( f, ph::_1, ph::_2, cref(x) );
}


} // namespace details 

using details::f_update_t;
using details::f_scalar_t;
using details::f_vector_t;
using details::f_dderiv_t;
using details::fix_direction;
using details::fix_position;


} // namespace ROL


#endif // ROL_FUNCTIONBINDING_HPP

