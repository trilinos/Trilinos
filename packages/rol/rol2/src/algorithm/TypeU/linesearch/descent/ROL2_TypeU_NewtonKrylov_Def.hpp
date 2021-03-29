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

#ifndef ROL2_TYPEU_NEWTONKRYLOV_DEF_HPP 
#define ROL2_TYPEU_NEWTONKRYLOV_DEF_HPP 

/** @ingroup step_group
    \class ROL2::TypeU::DescentDirection
    \brief Provides the interface to compute unconstrained optimization steps
           for line search.
*/

namespace ROL2 {
namespace TypeU {

template <class Real>
NewtonKrylov<Real>::HessianNK::HessianNK( const Ptr<Objective<Real>>&    obj,
                                         const Ptr<const Vector<Real>>& x ) 
  : obj_(obj), x_(x) {}

template <class Real>
void NewtonKrylov<Real>::HessianNK::apply(       Vector<Real>& Hv,
                                           const Vector<Real>& v,
                                                 Real&         tol ) const {
  obj_->hessVec(Hv,v,*x_,tol);
}

template <class Real>
NewtonKrylov<Real>::PrecondNK::PrecondNK( const Ptr<Objective<Real>>&    obj,
                                          const Ptr<const Vector<Real>>& x ) 
  : obj_(obj), x_(x) {
}

template <class Real>
void NewtonKrylov<Real>::PrecondNK::apply(       Vector<Real>& Hv,
                                           const Vector<Real>& v,
                                                 Real&         tol ) const {
  obj_->hessVec(Hv,v,*x_,tol);
}


template <class Real>
NewtonKrylov<Real>::NewtonKrylov( ParameterList& parlist ) {
  auto& klist = parlist.sublist("General").sublist("Krylov");
  auto& slist = parlist.sublist("General").sublist("Secant");
  useSecantPrecond_ = slist.get("Use as Preconditioner", false);
  krylovName_ = klist.get("Type","Conjugate Gradients");
  krylovType_ = Krylov<Real>::type_dict[krylovName_];
  krylov_ = Krylov<Real>::create(parlist);
  secantName_ = slist.get("Type","Limited-Memory BFGS");
  secantType_ = Secant<Real>::type_dict[secantName_];
  if ( useSecantPrecond_ ) {
    secant_  = Secant<Real>::create(parlist);
    precond_ = secant_;
  }
}
  
template <class Real>
NewtonKrylov<Real>::NewtonKrylov(       ParameterList&     parlist,
                                  const Ptr<Krylov<Real>>& krylov,
                                        Ptr<Secant<Real>>& secant,
                                       bool               computeObj )
  : secant_(secant), krylov_(krylov),
    krylovType_(Krylov<Real>::Type::UserDefined), 
    secantType_(Secant<Real>::Type::UserDefined) { 

  auto& klist = parlist.sublist("General").sublist("Krylov");
  auto& slist = parlist.sublist("General").sublist("Secant");

  // Initialize secant object
  if ( useSecantPrecond_ ) {
    if( secant_ != nullPtr ) {
      secantName_ = slist.get("Type","Limited-Memory BFGS");
      secantType_ = Secant<Real>::type_dict[secantName_];
    }
    else {
      secantName_ = slist.get("User Defined Secant Name",
                              "Unspecified User Defined Secant Method");
    }
    precond_ = secant_;
  }
  // Initialize Krylov object
  if ( krylov != nullPtr ) {
    krylovName_ = klist.get("Type","Conjugate Gradients");
    krylovType_ = Krylov<Real>::type_dict[krylovName_];
    krylov_ = Krylov<Real>::create(parlist);
  }
  else 
    krylovName_ = klist.get("User Defined Krylov Name",
                            "Unspecified User Defined Krylov Method");
}


template <class Real>
void NewtonKrylov<Real>::compute(       Vector<Real>&    s, 
                                        Real&            snorm, 
                                        Real&            sdotg, 
                                        int&             iter, 
                                        int&             flag,
                                  const Vector<Real>&    x, 
                                  const Vector<Real>&    g, 
                                        Objective<Real>& obj) {
  // Build Hessian and Preconditioner object
  auto obj_ptr = makePtrFromRef(obj);
  auto x_ptr   = makePtrFromRef(x);
  auto hessian = makePtr<HessianNK>(obj_ptr,x_ptr);
  Ptr<LinearOperator<Real>> precond;
  if ( !useSecantPrecond_ ) precond = makePtr<PrecondNK>(obj_ptr,x_ptr);

  // Run Krylov method
  flag = 0; iter = 0;
  krylov_->run(s,*hessian,g,*precond,iter,flag);

  // Check Krylov flags
  if ( flag == 2 && iter <= 1 )  s.set(g.dual());
  
  s.scale(static_cast<Real>(-1));
  snorm = s.norm();
  sdotg = s.apply(g);
}






template <class Real>
void NewtonKrylov<Real>::update( const Vector<Real>& x, 
                                const Vector<Real>& s,
                                const Vector<Real>& gold, 
                                const Vector<Real>& gnew,
                                      Real          snorm, 
                                      int           iter ) {
  if( useSecantPrecond_ ) 
    secant_->updateStorage(x,gnew,gold,s,snorm,iter+1);
}

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_DESCENTDIRECTION_H
