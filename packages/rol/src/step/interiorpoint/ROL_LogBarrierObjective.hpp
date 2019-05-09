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

#ifndef ROL_LOG_BARRIER_OBJECTIVE_H
#define ROL_LOG_BARRIER_OBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::LogBarrierObjective
    \brief Log barrier objective for interior point methods
*/

namespace ROL {

template <class Real>
class LogBarrierObjective : public Objective<Real> {
public:

    

  /* \brief Objective value J(x) = \f$-\sum_i \log(x_i) \f$ */
  Real value( const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<Vector<Real> > logx = x.clone();
    logx->set(x);
    
    Elementwise::Logarithm<Real> log;

    logx->applyUnary(log);

    Elementwise::ReductionSum<Real> sum;

    Real result = -(logx->reduce(sum));

    return result;
  }  

  /* \brief gradient g_i = \f$ 1/x_i \f$ */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    
    g.set(x);
     
    Elementwise::Reciprocal<Real> reciprocal;
   
    g.applyUnary(reciprocal);
    g.scale(-1.0);
  }

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
 
    ROL::Ptr<Vector<Real> > dbyx = d.clone();
    dbyx->set(x);

    struct Division : public Elementwise::BinaryFunction<Real> {
      Real apply( const Real &xc, const Real &dc ) const {
        return dc/xc;    
      }
    } division;

    dbyx->applyBinary( division, d );
    
    Elementwise::ReductionSum<Real> sum;

    return -dbyx->reduce(sum);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.set(v);
    
    struct HessianApply : public Elementwise::BinaryFunction<Real> {
      Real apply( const Real &vc, const Real &xc ) const { 
        return vc/(xc*xc);   
      }
    } hessian;

    hv.applyBinary(hessian,x);

  }

};

} // namespace ROL

#endif // ROL_LOG_BARRIER_OBJECTIVE_H
