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

#ifndef ROL_NEWTONKRYLOV_U_H
#define ROL_NEWTONKRYLOV_U_H

#include "ROL_DescentDirection_U.hpp"

#include "ROL_Types.hpp"
#include "ROL_Secant.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_LinearOperator.hpp"

/** @ingroup step_group
    \class ROL::NewtonKrylov_U
    \brief Provides the interface to compute optimization steps
           with projected inexact Newton's method using line search.
*/

namespace ROL {

template<typename Real>
class NewtonKrylov_U : public DescentDirection_U<Real> {
private:

  Ptr<Secant<Real>>         secant_; ///< Secant object (used for quasi-Newton)
  Ptr<Krylov<Real>>         krylov_; ///< Krylov solver object (used for inexact Newton)
  Ptr<LinearOperator<Real>> precond_;

  EKrylov ekv_;
  ESecant esec_;
 
  bool useSecantPrecond_; ///< Whether or not a secant approximation is used for preconditioning inexact Newton

  std::string krylovName_;
  std::string secantName_;

  class HessianNK : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>> obj_;
    const Ptr<const Vector<Real>> x_;
  public:
    HessianNK(const Ptr<Objective<Real>> &obj,
              const Ptr<const Vector<Real>> &x) : obj_(obj), x_(x) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      obj_->hessVec(Hv,v,*x_,tol);
    }
  };

  class PrecondNK : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>> obj_;
    const Ptr<const Vector<Real>> x_;
  public:
    PrecondNK(const Ptr<Objective<Real>> &obj,
              const Ptr<const Vector<Real>> &x) : obj_(obj), x_(x) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      obj_->precond(Hv,v,*x_,tol);
    }
  };

public:

  /** \brief Constructor.

      Standard constructor to build a NewtonKrylovStep object.  Algorithmic 
      specifications are passed in through a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  NewtonKrylov_U(ParameterList &parlist)
    : secant_(nullPtr), krylov_(nullPtr), useSecantPrecond_(false) {
    // Parse ParameterList
    ParameterList& Glist = parlist.sublist("General");
    useSecantPrecond_ = Glist.sublist("Secant").get("Use as Preconditioner", false);
    // Initialize Krylov object
    krylovName_ = Glist.sublist("Krylov").get("Type","Conjugate Gradients");
    ekv_ = StringToEKrylov(krylovName_);
    krylov_ = KrylovFactory<Real>(parlist);
    // Initialize secant object
    secantName_ = Glist.sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_ = StringToESecant(secantName_);
    if ( useSecantPrecond_ ) {
      secant_  = SecantFactory<Real>(parlist);
      precond_ = secant_;
    }
  }

  /** \brief Constructor.

      Constructor to build a NewtonKrylovStep object with user-defined 
      secant and Krylov objects.  Algorithmic specifications are passed in through 
      a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     krylov     is a user-defined Krylov object
      @param[in]     secant     is a user-defined secant object
  */
  NewtonKrylov_U(ParameterList &parlist, const Ptr<Krylov<Real>> &krylov,
                 Ptr<Secant<Real>> &secant, const bool computeObj = true)
    : secant_(secant), krylov_(krylov),
      ekv_(KRYLOV_USERDEFINED), esec_(SECANT_USERDEFINED),
      useSecantPrecond_(false) {
    // Parse ParameterList
    ParameterList& Glist = parlist.sublist("General");
    useSecantPrecond_ = Glist.sublist("Secant").get("Use as Preconditioner", false);
    // Initialize secant object
    if ( useSecantPrecond_ ) {
      if(secant_ == nullPtr ) {
        secantName_ = Glist.sublist("Secant").get("Type","Limited-Memory BFGS");
        esec_ = StringToESecant(secantName_);
        secant_ = SecantFactory<Real>(parlist);
      }
      else {
        secantName_ = Glist.sublist("Secant").get("User Defined Secant Name",
                                                  "Unspecified User Defined Secant Method");
      }
      precond_ = secant_;
    }
    // Initialize Krylov object
    if ( krylov_ == nullPtr ) {
      krylovName_ = Glist.sublist("Krylov").get("Type","Conjugate Gradients");
      ekv_ = StringToEKrylov(krylovName_);
      krylov_ = KrylovFactory<Real>(parlist);
    }
    else {
      krylovName_ =  Glist.sublist("Krylov").get("User Defined Krylov Name",
                                                 "Unspecified User Defined Krylov Method");
    }
  }

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    // Build Hessian and Preconditioner object
    Ptr<Objective<Real>>    obj_ptr = makePtrFromRef(obj);
    Ptr<const Vector<Real>> x_ptr   = makePtrFromRef(x);
    Ptr<LinearOperator<Real>> hessian
      = makePtr<HessianNK>(obj_ptr,x_ptr);
    Ptr<LinearOperator<Real>> precond;
    if ( !useSecantPrecond_ ) {
      precond = makePtr<PrecondNK>(obj_ptr,x_ptr);
    }

    // Run Krylov method
    flag = 0; iter = 0;
    krylov_->run(s,*hessian,g,*precond,iter,flag);

    // Check Krylov flags
    if ( flag == 2 && iter <= 1 ) {
      s.set(g.dual());
    }
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    //sdotg = s.dot(g.dual());
    sdotg = s.apply(g);
  }

  void update(const Vector<Real> &x, const Vector<Real> &s,
              const Vector<Real> &gold, const Vector<Real> &gnew,
              const Real snorm, const int iter) override {
    // Update Secant Information
    if ( useSecantPrecond_ ) {
      secant_->updateStorage(x,gnew,gold,s,snorm,iter+1);
    }
  }

  std::string printName(void) const override {
    std::stringstream name;
    name << "Newton-Krylov Method using " << krylovName_;
    if (useSecantPrecond_) {
      name << " with " << secantName_ << " preconditioning";
    }
    return name.str();
  }
}; // class NewtonKrylov_U

} // namespace ROL

#endif
