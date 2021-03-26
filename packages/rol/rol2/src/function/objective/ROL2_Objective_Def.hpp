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

#ifndef ROL2_OBJECTIVE_DEF_H
#define ROL2_OBJECTIVE_DEF_H


namespace ROL2 {


template<typename Real>
Real Objective<Real>::dirDeriv( const Vector<Real>& x, 
                                const Vector<Real>& d, 
                                      Real&         tol ) {

  Real one(1);
  auto f = [this,&tol]( const Vector<Real>& x ) { return this->value(x,tol); };
  auto u = [this]( const Vector<Real>& x ) { this->update(x); };
  Real h = std::max(one,x.norm()/d.norm())*tol;
  return UFD::diff_scalar(f,u,x,d,h,1);
} // dirDeriv


template<typename Real>
void Objective<Real>::gradient(       Vector<Real>& g, 
                                const Vector<Real>& x, 
                                      Real&         tol ) {
  g.zero();
  
  Real deriv = 0, h = 0, xi = 0;

  for( int i=0; i<x.dimension(); ++i) {
    xi    = std::abs(x.dot(*x.basis(i)));
    h     =  ((xi < ROL_EPSILON<Real>) ? 1. : xi)*tol;
    deriv = this->dirDeriv(x,*x.basis(i),h);
    g.axpy(deriv,*g.basis(i));   
  }
} // gradient


template<typename Real>
void Objective<Real>::hessVec(       Vector<Real>& hv, 
                               const Vector<Real>&   v, 
                               const Vector<Real>&   x, 
                                     Real&           tol ) {
  
  Real zero(0), one(1);
  if ( v.norm() == zero )  hv.zero();
  else {

    Real h = std::max(one,x.norm()/v.norm())*tol;
    auto update = [this]( const auto& x ) { this->update(x); };
    auto fv     = [this,&tol] ( auto& g, const auto& x ) { this->gradient(g,x,tol); }; 
    UFD::diff_vector(fv,update,hv,x,v,h,1);
  }
} 



template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkGradient( const Vector<Real>& x,
                                const Vector<Real>& g,
                                const Vector<Real>& d,
                                      bool          printToStream,
                                      std::ostream& outStream,
                                      int           numSteps,
                                      int           order ) {
  Real tol = default_tolerance<Real>();

  auto f_value    = [this,&tol] ( const auto& x ) { return this->value(x,tol);  };
  auto f_gradient = [this,&tol] ( auto& g, const auto& x ) {  this->gradient(g,x,tol); };
  auto f_update   = [this] ( const auto& x ) {  this->update(x);  };
  
  ValidateFunction<Real> validator(order,numSteps,11,printToStream,outStream);
 
  return validator.scalar_derivative_check(f_value,f_gradient,f_update,g,d,x,"grad'*dir"); 
} // checkGradient



template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkGradient( const Vector<Real>&      x,
                                const Vector<Real>&      g,
                                const Vector<Real>&      d,
                                const std::vector<Real>& steps,
                                      bool               printToStream,
                                      std::ostream&      outStream,
                                      int                order ) {
 
  Real tol = default_tolerance<Real>();

  auto f_update   = [this] ( const auto& x ) {  this->update(x); };
  auto f_value    = [this,&tol] ( const auto& x ) { return this->value(x,tol); };
  auto f_gradient = [this,&tol] ( auto& g, const auto& x ) { this->gradient(g,x,tol);  };

  
  ValidateFunction<Real> validator(steps,order,11,printToStream,outStream);
 
  return validator.scalar_derivative_check(f_value,f_gradient,f_update,g,d,x,"grad'*dir");
} // checkGradient



template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkGradient( const Vector<Real>& x,
                                const Vector<Real>& d,
                                      bool          printToStream,
                                      std::ostream& outStream,
                                      int           numSteps,
                                      int           order ) { 
  return checkGradient(x, x.dual(), d, printToStream, outStream, numSteps, order);
}


template<typename Real>
std::vector<std::vector<Real>>
Objective<Real>::checkGradient( const Vector<Real>&      x,
                                const Vector<Real>&      d,
                                const std::vector<Real>& steps,
                                      bool               printToStream,
                                      std::ostream&      outStream,
                                      int                order ) {
  return checkGradient(x, x.dual(), d, steps, printToStream, outStream, order);
}







template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkHessVec( const Vector<Real>& x,
                               const Vector<Real>& hv,
                               const Vector<Real>& v,
                                     bool          printToStream,
                                     std::ostream& outStream,
                                     int           numSteps,
                                     int           order ) {
  std::vector<Real> steps(numSteps);
  for(int i=0;i<numSteps;++i)  steps[i] = pow(10,-i);
  return checkHessVec(x,hv,v,steps,printToStream,outStream,order);
} // checkHessVec



template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkHessVec( const Vector<Real>&      x,
                               const Vector<Real>&      hv,
                               const Vector<Real>&      v,
                               const std::vector<Real>& steps,
                                     bool               printToStream,
                                     std::ostream&      outStream,
                                     int                order ) {

  Real tol = default_tolerance<Real>();

  auto f_gradient = [this,&tol] ( auto& g, const auto& x ) { this->gradient(g,x,tol);};
  auto f_hessVec  = [this,&tol] ( auto& hv, const auto& v, const auto& x ) { this->hessVec(hv,v,x,tol); };
  auto f_update   = [this] ( const Vector<Real>& x ) { this->update(x); };
  
  ValidateFunction<Real> validator(steps,order,11,printToStream,outStream);
 
  return validator.vector_derivative_check(f_gradient,f_hessVec,f_update,hv,v,x,"norm(Hess*vec)");

} // checkHessVec


template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkHessVec( const Vector<Real>& x,
                               const Vector<Real>& v,
                                     bool          printToStream,
                                     std::ostream& outStream,
                                     int           numSteps,
                                     int           order ) {

  return checkHessVec(x, x.dual(), v, printToStream, outStream, numSteps, order);
}


template<typename Real>
std::vector<std::vector<Real>> 
Objective<Real>::checkHessVec( const Vector<Real>&      x,
                               const Vector<Real>&      v,
                               const std::vector<Real>& steps,
                                     bool               printToStream,
                                     std::ostream&      outStream,
                                     int                order ) {
  return checkHessVec(x, x.dual(), v, steps, printToStream, outStream, order);
}

template<typename Real>
std::vector<Real> 
Objective<Real>::checkHessSym( const Vector<Real>& x,
                               const Vector<Real>& hv,
                               const Vector<Real>& v,
                               const Vector<Real>& w,
                                     bool          printToStream,
                                     std::ostream& outStream ) {

  Real tol = default_tolerance<Real>();
  
  // Compute (Hessian at x) times (vector v).
  auto h = hv.clone();
  this->hessVec(*h, v, x, tol);
  Real wHv = w.dot(h->dual());

  this->hessVec(*h, w, x, tol);
  Real vHw = v.dot(h->dual());

  std::vector<Real> hsymCheck(3, 0);

  hsymCheck[0] = wHv;
  hsymCheck[1] = vHw;
  hsymCheck[2] = std::abs(vHw-wHv);

  // Save the format state of the original outStream.
  NullStream oldFormatState;
  oldFormatState.copyfmt(outStream);

  if (printToStream) {
    outStream << std::right
              << std::setw(20) << "<w, H(x)v>"
              << std::setw(20) << "<v, H(x)w>"
              << std::setw(20) << "abs error"
              << "\n";
    outStream << std::scientific << std::setprecision(11) << std::right
              << std::setw(20) << hsymCheck[0]
              << std::setw(20) << hsymCheck[1]
              << std::setw(20) << hsymCheck[2]
              << "\n";
  }

  // Reset format state of outStream.
  outStream.copyfmt(oldFormatState);

  return hsymCheck;

} // checkHessSym

template<typename Real>
std::vector<Real> 
Objective<Real>::checkHessSym( const Vector<Real>& x,
                               const Vector<Real>& v,
                               const Vector<Real>& w,
                                     bool          printToStream,
                                     std::ostream& outStream ) {
  return checkHessSym(x, x.dual(), v, w, printToStream, outStream);
}


} // namespace ROL2

#endif // ROL2_OBJECTIVE_DEF_HPP

