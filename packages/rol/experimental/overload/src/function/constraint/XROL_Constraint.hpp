
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

#include "XROL_ConstraintVectors.hpp"


namespace XROL {


template<class X, class C>
class Constraint {

public:

  Constraint();

  virtual ~Constraint();

  virtual void update( const X& x );

  virtual void value( C& c, const X& x, magnitude_t<X>& tol ) = 0;

  virtual void applyJacobian( C& jv, 
                              const X& v, 
                              const X& x, 
                              magnitude_t<X>& tol );

  virtual void applyAdjointJacobian( dual_t<X>& ajv, 
                                     const dual_t<C>& v,
                                     const X& x,
                                     const C& vdual,
                                     magnitude_t<X>& tol );


  virtual void applyAdjointHessian( dual_t<X>& ahuv,
                                    const dual_t<C>& u,
                                    const X& v,
                                    const X& x,
                                    magnitude_t<X>& tol );   
/*
  virtual std::vector<magnitude_t<X>>
  solveAugmentedSystem( X&               v1,
                        dual_t<C>&       v2,
                        const dual_t<C>& b1,
                        const C&         b2,
                        const X&         x,
                        magnitude_t<X>&  tol );
*/

  virtual void applyPreconditioner( dual_t<C>& pv,
                                    const C& v, 
                                    const X& x,
                                    const dual_t<X>& g,
                                    magnitude_t<X>& tol );
  virtual void activate( void );
 
  virtual void deactivate( void );

private:

  bool activated_;

  // Constraint space scratch vector in the event of default 
  // adjoint Jacobian use
  std::unique_ptr<ConstraintVectors<X,C>> convec_; 
  
};


} // namespace XROL

