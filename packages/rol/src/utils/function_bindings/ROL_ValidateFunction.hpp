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
#ifndef ROL_VALIDATEFUNCTION_HPP
#define ROL_VALIDATEFUNCTION_HPP

/** \file  ROL_ValidateFunction.hpp
    \brief Provides a set of tools for validating the behavior of several
           function types that are commonly used in ROL.
 
           - Finite difference check of derivatives
           - Symmetry check of linear operators
           - Adjoint consistency check for linear operators
           - Inverse identity check for linear operators
*/

#include "ROL_FiniteDifference.hpp"

namespace ROL {

namespace details {

template<typename Real>
class ValidateFunction {
public:

  using V = Vector<Real>;

  ValidateFunction( const int order = 1,
                    const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                    const int width = 20,
                    const int precision = 11,
                    const bool printToStream = true,
                    ostream& os = cout );  

  virtual ~ValidateFunction(){}

  virtual vector<vector<Real>> derivative_check( f_scalar_t<Real> f_value, 
                                                 f_vector_t<Real> f_derivative, 
                                                 f_update_t<Real> f_update,
                                                 const V& g, 
                                                 const V& x,
                                                 const V& v,
                                                 const string& label ) const;
                                        
  virtual vector<vector<Real>> derivative_check( f_vector_t<Real> f_value, 
                                                 f_dderiv_t<Real> f_derivative, 
                                                 f_update_t<Real> f_update,
                                                 const V& c, 
                                                 const V& x,
                                                 const V& v,
                                                 const string& label ) const;

  virtual vector<Real> symmetry_check( f_dderiv_t<Real> A, 
                                       f_update_t<Real> A_update,
                                       const V& u, 
                                       const V& v, 
                                       const V& x,
                                       const string& name="Linear Operator",
                                       const string& symbol="A" ) const;

  virtual vector<Real> adjoint_consistency_check( f_dderiv_t<Real> A,
                                                  f_dderiv_t<Real> A_adj,
                                                  f_update_t<Real> A_update,
                                                  const V& u,  
                                                  const V& v,
                                                  const V& x, 
                                                  const string& name="Linear Operator",
                                                  const string& symbol="A" ) const;

  virtual vector<Real> inverse_check( f_dderiv_t<Real> A,
                                      f_dderiv_t<Real> A_inv,
                                      f_update_t<Real> A_update,
                                      const V& v,
                                      const V& x, 
                                      const string& name="Linear Operator",
                                      const string& symbol="A" ) const;
private:

  int                    order_;         // Finite difference order
  int                    numSteps_;      // Number of evalutions of different step sizes
  int                    width_;         // For print formatting
  int                    precision_;     // Number of digits to display
  bool                   printToStream_; // False will suppress output
  vector<Real>           steps_;         // Set of step sizes of FD approximation
  
  ostream& os_;                          // pointer to Output stream
 
  mutable Ptr<VectorWorkspace<Real>> workspace_;
  FiniteDifference<Real>             fd_;


}; // ValidateFunction


} // namespace details 

using details::ValidateFunction;

} // namespace ROL

#include "ROL_ValidateFunctionDef.hpp"

#endif // ROL_CHECKFUNCTION_HPP

