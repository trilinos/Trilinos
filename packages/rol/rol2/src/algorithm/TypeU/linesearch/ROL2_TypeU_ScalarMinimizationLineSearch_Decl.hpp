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

#ifndef ROL2_TYPEU_SCALARMINIMIZATIONLINESEARCH_DECL_H
#define ROL2_TYPEU_SCALARMINIMIZATIONLINESEARCH_DECL_H

/** \class ROL2::TypeU::ScalarMinimizationLineSearch
    \brief Implements line search methods that attempt to minimize the
           scalar function \f$\phi(t) := f(x+ts)\f$.
*/

namespace ROL2 {
namespace TypeU {

template<typename Real>
class ScalarMinimizationLineSearch : public LineSearch<Real> {
public:

  //----------------------------------------------------------

  class Phi : public ScalarFunction<Real> {
  public:
    Phi( const Ptr<Vector<Real>>&       xnew,
         const Ptr<const Vector<Real>>& x,
         const Ptr<const Vector<Real>>& s,
         const Ptr<Objective<Real>>&    obj,
               bool                     FDdirDeriv = false );

  Real value( Real alpha ) override;

  Real deriv(const Real alpha) override;

  private:
    const Ptr<Vector<Real>>       xnew_;
    const Ptr<const Vector<Real>> x_, s_;
    const Ptr<Objective<Real>>    obj_;
    Real ftol_, alpha_, val_;
    bool FDdirDeriv_;
  }; // class Phi
 
  //----------------------------------------------------------

  class StatusTest : public ScalarMinimizationStatusTest<Real> {
  public:
    StatusTest( Real                             f0, 
                Real                             g0,
                Real                             c1, 
                Real                             c2, 
                Real                             c3,
                int                              max_nfval, 
                LineSearch<Real>::CurvartureCond ccond,
                const Ptr<ScalarFunction<Real>>& phi)
      : phi_(phi), f0_(f0), g0_(g0), c1_(c1), c2_(c2), c3_(c3),
        max_nfval_(max_nfval), ccond_(econd) {}

    bool check( Real& x, 
                Real& fx, 
                Real& gx,
                int&  nfval, 
                int&  ngval, 
                bool  deriv = false ) override;
  private:
    Ptr<ScalarFunction<Real>> phi_;
    const Real f0_, g0_, c1_, c2_, c3_;
    const int max_nfval_;
    const LineSearch<Real>::CurvatureConditionU ccond_;
  }; // class StatusTest

  //----------------------------------------------------------

  // Constructor
  ScalarMinimizationLineSearch(       ParameterList&                 parlist, 
                                const Ptr<ScalarMinimization<Real>>& sm = nullPtr,
                                const Ptr<Bracketing<Real>>&         br = nullPtr,
                                const Ptr<ScalarFunction<Real>>&     sf = nullPtr );
    
  void initialize( const Vector<Real>& x,
                   const Vector<Real>& g ) override {
    LineSearch<Real>::initialize(x,g);
    xnew_ = x.clone();
    g_    = g.clone();
  }

  // Find the minimum of phi(alpha) = f(x + alpha*s) using Brent's method
  void run(       Real&            alpha,
                  Real&            fval,
                  int&             ls_neval,
                  int&             ls_ngrad,
            const Real&            gs,
            const Vector<Real>&    s,
            const Vector<Real>&    x,
                  Objective<Real>& obj ) override;
private:

  Ptr<Vector<Real>>             xnew_; 
  Ptr<Vector<Real>>             g_;
  Ptr<ScalarMinimization<Real>> sm_;
  Ptr<Bracketing<Real>>         br_;
  Ptr<ScalarFunction<Real>>     sf_;

  Type ccond_;
  Real c1_;
  Real c2_;
  Real c3_;
  int  max_nfval_;
  bool FDdirDeriv_;


}; 

} // namespace TypeU
} // namespace ROL

#endif // ROL2_TYPEU_SCALARMINIMIZATIONLINESEARCH_DECL_H
