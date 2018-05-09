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
#ifndef ROL_OBJECTIVE_CHECKINTERFACE_HPP
#define ROL_OBJECTIVE_CHECKINTERFACE_HPP

#include "ROL_Objective.hpp"
#include "ROL_FunctionBindings.hpp"

namespace ROL {
namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class Objective_CheckInterface {
private:
  using V = Vector<Real>;
  Objective<Real>& obj_;
  Real tol_;

public:

  Objective_CheckInterface( Objective<Real>& obj ) : 
    obj_(obj), tol_(sqrt(ROL_EPSILON<Real>())) {}
   
  f_update_t<Real> update() {
    return bind( &Objective<Real>::update, &obj_, ph::_1, true, 0 );
  }

  f_scalar_t<Real> value() {
    return bind( &Objective<Real>::value, &obj_, ph::_1, tol_);
  }

  f_vector_t<Real> gradient() {
    return bind( &Objective<Real>::gradient, &obj_, ph::_1, ph::_2, tol_);
  }

  f_dderiv_t<Real> hessVec() {
    return bind( &Objective<Real>::hessVec, &obj_, ph::_1, ph::_2, ph::_3, tol_);
  }

}; // Objective_CheckInterface

} // namespace details

using details::Objective_CheckInterface;
template<typename Real>
Objective_CheckInterface<Real> make_check( Objective<Real>& obj ) {
  return Objective_CheckInterface<Real>(obj);
}

} // namespace ROL


#endif // ROL_OBJECTIVE_CHECKINTERFACE_HPP

