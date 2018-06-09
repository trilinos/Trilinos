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
#ifndef ROL_FINITEDIFFERENCE_HPP
#define ROL_FINITEDIFFERENCE_HPP


/** \class ROL::FiniteDifference
    \brief General interface for numerical approximation of Constraint 
           or Objective methods using the finite difference method. 

           TODO: Should optionally store values for reuse as needed
*/

#include "ROL_FunctionBindings.hpp"

#include "ROL_Vector.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace ROL {

namespace details {

using namespace std;

template<typename Real>
class FiniteDifference {
public:
 
  using V = ROL::Vector<Real>;

  FiniteDifference( const int order = 1 );
  FiniteDifference( const int order, const Ptr<VectorWorkspace<Real>>& workspace );

  virtual ~FiniteDifference(){}

  /** Approximately compute the derivative of f(x) in the direction v 
      using the step size h */
  virtual Real operator()( f_scalar_t<Real>& f_value,
                           f_update_t<Real>& f_update,
                           const V& v, 
                           const V& x,
                           const Real h ) const;

  /** Approximately compute the Jacobian of f(x) applied to the direction v
      using the step size h */
  virtual void operator()( f_vector_t<Real>& f_value, 
                           f_update_t<Real>& f_update,
                           V& Jv, 
                           const V& v, 
                           const V& x,
                           const Real h ) const;             
private:

  const int order_; // Finite difference order (1,2,3, or 4)
  
  mutable Ptr<VectorWorkspace<Real>> workspace_;

}; // class FiniteDifference

} // namespace details

using details::FiniteDifference;


} // namespace ROL

#include "ROL_FiniteDifferenceDef.hpp"

#endif // ROL_FINITEDIFFERENCE_HPP

