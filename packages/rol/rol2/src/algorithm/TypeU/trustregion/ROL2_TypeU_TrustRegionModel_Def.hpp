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

#ifndef ROL2_TYPEU_TRUSTREGIONMODEL_DEF_H
#define ROL2_TYPEU_TRUSTREGIONMODEL_DEF_H

namespace ROL2 {
namespace TypeU {


template<typename Real>
TrustRegionModel<Real>::TrustRegionModel(       ParameterList&     list,
                                          const Ptr<Secant<Real>>& secant,
                                                Secant<Real>::Mode mode ) : 
    : obj_(nullPtr), x_(nullPtr), g_(nullPtr), secant_(secant) {
    auto& slist = list.sublist("General").sublist("Secant");
    useSecantPrecond_ = slist.get("Use as Preconditioner", false);
    useSecantHessVec_ = slist.get("Use as Hessian",        false);
    if (secant_ == nullPtr) secant_ = SecantFactory<Real>(list,mode);
  }

template<typename Real>
void TrustRegionModel<Real>::initialize( const Vector<Real>& x, 
                   const Vector<Real>& g ) {
    dual_ = g.clone();
  }

  // Some versions of Clang will issue a warning that update hides and 
  // overloaded function without this using declaration
  using Objective<Real>::update;

template<typename Real>
void TrustRegionModel<Real>::validate( Objective<Real>&    obj,
                 const Vector<Real>& x,
                 const Vector<Real>& g,
                ETrustRegionU       etr) {
    if ( !useSecantHessVec_ &&
        (etr == TrustRegion<Real>::Type::DogLeg
       || etr == TrustRegion<Real>::Type::DoubleDogLeg) {
      try {
        Real htol = default_tolerance<Real>();
        auto v  = g.clone();
        auto hv = x.clone();
        obj.invHessVec(*hv,*v,x,htol);
      }
      catch (std::exception &e) {
        useSecantHessVec_ = true;
      }
    }
  }

template<typename Real>
void TrustRegionModel<Real>::setData(       Objective<Real>& obj,
                                      const Vector<Real>&    x,
                                      const Vector<Real>&    g) {
  obj_ = makePtrFromRef(obj);
  x_   = makePtrFromRef(x);
  g_   = makePtrFromRef(g);
}

template<typename Real>
void TrustRegionModel<Real>::update( const Vector<Real>& x, 
                                     const Vector<Real>& s,
                                     const Vector<Real>& gold, 
                                     const Vector<Real>& gnew,
                                           Real          snorm, 
                                           int           iter) {
  // Update Secant Information
  if (useSecantHessVec_ || useSecantPrecond_) {
    secant_->updateStorage(x,gnew,gold,s,snorm,iter);
  }
}

template<typename Real>
Real TrustRegionModel<Real>::value( const Vector<Real> &s, Real &tol ) {
  applyHessian(*dual_,s,tol);
  dual_->scale(static_cast<Real>(0.5));
  dual_->plus(*g_);
  return dual_->apply(s);
}

template<typename Real>
void TrustRegionModel<Real>::gradient(       Vector<Real>& g, 
                                       const Vector<Real>& s, 
                                             Real&         tol ) {
  applyHessian(g,s,tol);
  g.plus(*g_);
}

template<typename Real>
void TrustRegionModel<Real>::applyHessian(       Vector<Real>& hv,
                                           const Vector<Real>& v, 
                                                 Real&         tol) {
  if ( useSecantHessVec_ && secant_ != nullPtr )  secant_->applyB(hv,v);
  else obj_->hessVec(hv,v,*x_,tol);
}

template<typename Real>
void TrustRegionModel<Real>::applyInvHessian(       Vector<Real>& hv, 
                                              const Vector<Real>& v, 
                                                    Real&         tol) {
  if ( useSecantHessVec_ && secant_ != nullPtr ) secant_->applyH(hv,v);
  else obj_->invHessVec(hv,v,*x_,tol);
}

template<typename Real>
void TrustRegionModel<Real>::applyPrecond(       Vector<Real>& Pv, 
                                           const Vector<Real>& v, 
                                                 Real&         tol ) {
  if ( useSecantPrecond_  && secant_ != nullPtr )  secant_->applyH(Pv,v);
  else  obj_->precond(Pv,v,*x_,tol);
}

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_TRUSTREGIONMODEL_DEF_H
