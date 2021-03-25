// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL2) Package
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
#ifndef ROL2_LSR1_DECL_HPP
#define ROL2_LSR1_DECL_HPP

/** \class ROL2::lSR1
    \brief Provides interface for and implements limited-memory SR1 operators.
*/

namespace ROL2 {

template<class Real>
class lSR1 : public LinearOperator<Real> {
public:

  ~lSR1() = default;

  // Constructor
  lSR1( int                M = 10, 
        bool               useDefaultScaling = true, 
        Real               Bscaling = Real(1), 
        Secant<Real>::Mode mode = Secant<Real>::Mode::Both ) {
    if( useDefaultScaling_ ) Bscaling_ = 1;
  }

  // Update lSR1 Approximation
  void updateStorage( const Vector<Real>& x,  
                      const Vector<Real>& grad,
                      const Vector<Real>& gp, 
                      const Vector<Real>& s,
                            Real          snorm,      
                            int           iter ) override;

  // Apply lSR1 Approximate Inverse Hessian
  void applyH(       Vector<Real>& Hv, 
               const Vector<Real>& v ) const override;

  // Apply Initial lSR1 Approximate Inverse Hessian
  void applyH0(       Vector<Real>& Hv, 
                const Vector<Real>& v ) const override;

  // Apply lSR1 Approximate Hessian
  void applyB(       Vector<Real>& Bv, 
               const Vector<Real>& v ) const override;

  // Apply Initial lSR1 Approximate Hessian 
  void applyB0(       Vector<Real>& Bv, 
                const Vector<Real>& v ) const override;

private:

  Ptr<Vector<Real>> Bs_, Hy_, prim_, dual_;
  mutable bool H0called_, B0called_;
  bool isInitialized_;

  using Secant<Real>::y_;
  using Secant<Real>::useDefaultScaling_;
  using Secant<Real>::Bscaling_;

}; // lSR1

} // namespace ROL2

#endif // ROL2_LSR1_DECL_HPP
