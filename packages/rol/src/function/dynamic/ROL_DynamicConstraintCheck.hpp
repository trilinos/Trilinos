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
#ifndef ROL_DYNAMICCONSTRAINTCHECK_HPP
#define ROL_DYNAMICCONSTRAINTCHECK_HPP

#include <functional>

#include "ROL_DynamicConstraint.hpp"

namespace ROL {

namespace details {

using namespace std;
using Finite_Difference_Arrays::shifts;
using Finite_Difference_Arrays::weights;

template<typename Real>
class DynamicConstraintCheck : public DynamicConstraint<Real> {
public:

  using V  = Vector<Real>;  
  using DC = DynamicConstraint<Real>;
 

  // Provide vectors from which to create random direction vectors
  // (No bound constraint)
  DynamicConstraintCheck( DynamicConstraint<Real>& con,
                          ROL::ParameterList& pl, 
                          ostream& os );

  DynamicConstraintCheck( DynamicConstraint<Real>& con,
                          const int numSteps, const int order,
                          ostream& os );
 

  virtual ~DyanamicConstraintCheck() {}

private:

  DC<Real>&    con_;
  ostream&     os_;   

  int          order_;
  int          numSteps_;
  vector<Real> steps_;
};

} // namespace details

using details::DynamicConstraintCheck;

} // namespace ROL

#include "ROL_DynamicConstraintCheckDef.hpp"

#endif // ROL_DYNAMICCONSTRAINTCHECK_HPP


