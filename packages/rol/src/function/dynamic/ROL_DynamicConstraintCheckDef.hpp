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
#ifndef ROL_DYNAMICCONSTRAINTCHECKDEF_HPP
#define ROL_DYNAMICCONSTRAINTCHECKDEF_HPP


namespace ROL {

namespace details {

using namespace std;
namespace ph = std::placeholders; // defines ph::_1, ph::_2, ...

template<typename Real>
DynamicConstraintCheck<Real>::DynamicConstraintCheck( DynamicConstraint<Real>& con,
                                                      ROL::ParameterList& pl,
                                                      ostream& os=cout ) :
  con_(con), os_(os) {
  auto& fdlist = pl.sublist("Finite Difference Check"); 
  order_ = fdlist.get("Order", 1);
  numSteps_ = fdlist.get("Number of Steps", 
}

template<typename Real>
DynamicConstraintCheck<Real>::DynamicConstraintCheck( DynamicConstraint<Real>& con,
                                                     const int numSteps = ROL_NUM_CHECKDERIVSTEPS, 
                                                     const int order = 1,
                                                     ostream& os=cout ) : 
  con_(con), numSteps_(numSteps), order_(order), os_(os) {
  steps_.resize(numSteps_);
}



template<typename Real>
DynamicConstraintCheck<Real>::value( V& c, const V& uo, const V& un, 
                                     const V& z, const TS& ts ) const override {
  con_.value(c,uo,un,z,ts);
}
 
template<typename Real>
DynamicConstraintCheck<Real>::applyJacobian_uo( V& jv, const V& vo, const V& uo, 
                                                const V& un, const V& z, 
                                                const TS& ts ) const override {

  auto f_val = bind( &DC::value, ph::_2, ph::_1, un, z, ts );
  auto f_der = bind( &DC::applyJacobian_uo, ph::_3, ph::_2, ph::_1, un, z, ts );
   
}


} // namespace details

} // namespace ROL


#endif // ROL_DYNAMICCONSTRAINTCHECKDEF_HPP

