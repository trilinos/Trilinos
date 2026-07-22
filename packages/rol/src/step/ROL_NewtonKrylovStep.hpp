// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NEWTONKRYLOVSTEP_H
#define ROL_NEWTONKRYLOVSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"

#include "ROL_Secant.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_LinearOperator.hpp"

#include <sstream>
#include <iomanip>

/** @ingroup step_group
    \class ROL::NewtonKrylovStep
    \brief Provides the interface to compute optimization steps
           with projected inexact Newton's method using line search.
*/

namespace ROL {

template <class Real>
class NewtonKrylovStep : public Step<Real> {
private:

  ROL::Ptr<Secant<Real> > secant_; ///< Secant object (used for quasi-Newton)
  ROL::Ptr<Krylov<Real> > krylov_; ///< Krylov solver object (used for inexact Newton)

  EKrylov ekv_;
  ESecant esec_;

  ROL::Ptr<Vector<Real> > gp_;
 
  int iterKrylov_; ///< Number of Krylov iterations (used for inexact Newton)
  int flagKrylov_; ///< Termination flag for Krylov method (used for inexact Newton)
  int verbosity_;  ///< Verbosity level
  const bool computeObj_;
 
  bool useSecantPrecond_; ///< Whether or not a secant approximation is used for preconditioning inexact Newton

  std::string krylovName_;
  std::string secantName_;


  class HessianNK : public LinearOperator<Real> {
  private:
    const ROL::Ptr<Objective<Real> > obj_;
    const ROL::Ptr<Vector<Real> > x_;
  public:
    HessianNK(const ROL::Ptr<Objective<Real> > &obj,
              const ROL::Ptr<Vector<Real> > &x) : obj_(obj), x_(x) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      obj_->hessVec(Hv,v,*x_,tol);
    } 
  };

  class PrecondNK : public LinearOperator<Real> {
  private:
    const ROL::Ptr<Objective<Real> > obj_;
    const ROL::Ptr<Vector<Real> > x_;
  public:
    PrecondNK(const ROL::Ptr<Objective<Real> > &obj,
              const ROL::Ptr<Vector<Real> > &x) : obj_(obj), x_(x) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      obj_->precond(Hv,v,*x_,tol);
    }
  };

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  /** \brief Constructor.

      Standard constructor to build a NewtonKrylovStep object.  Algorithmic 
      specifications are passed in through a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  NewtonKrylovStep( ROL::ParameterList &parlist, const bool computeObj = true )
    : Step<Real>(), secant_(ROL::nullPtr), krylov_(ROL::nullPtr),
      gp_(ROL::nullPtr), iterKrylov_(0), flagKrylov_(0),
      verbosity_(0), computeObj_(computeObj), useSecantPrecond_(false) {
    // Parse ParameterList
    ROL::ParameterList& Glist = parlist.sublist("General");
    useSecantPrecond_ = Glist.sublist("Secant").get("Use as Preconditioner", false);
    verbosity_ = Glist.get("Print Verbosity",0);
    // Initialize Krylov object
    krylovName_ = Glist.sublist("Krylov").get("Type","Conjugate Gradients");
    ekv_ = StringToEKrylov(krylovName_);
    krylov_ = KrylovFactory<Real>(parlist);
    // Initialize secant object
    secantName_ = Glist.sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_ = StringToESecant(secantName_);
    if ( useSecantPrecond_ ) {
      secant_ = SecantFactory<Real>(parlist);
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
  NewtonKrylovStep(ROL::ParameterList &parlist,
             const ROL::Ptr<Krylov<Real> > &krylov,
             const ROL::Ptr<Secant<Real> > &secant,
             const bool computeObj = true)
    : Step<Real>(), secant_(secant), krylov_(krylov),
      ekv_(KRYLOV_USERDEFINED), esec_(SECANT_USERDEFINED),
      gp_(ROL::nullPtr), iterKrylov_(0), flagKrylov_(0),
      verbosity_(0), computeObj_(computeObj), useSecantPrecond_(false) {
    // Parse ParameterList
    ROL::ParameterList& Glist = parlist.sublist("General");
    useSecantPrecond_ = Glist.sublist("Secant").get("Use as Preconditioner", false);
    verbosity_ = Glist.get("Print Verbosity",0);
    // Initialize secant object
    if ( useSecantPrecond_ ) {
      if(secant_ == ROL::nullPtr ) {
        secantName_ = Glist.sublist("Secant").get("Type","Limited-Memory BFGS");
        esec_ = StringToESecant(secantName_);
        secant_ = SecantFactory<Real>(parlist);
      }
      else {
        secantName_ = Glist.sublist("Secant").get("User Defined Secant Name",
                                                  "Unspecified User Defined Secant Method");         
      }
    }
    // Initialize Krylov object
    if ( krylov_ == ROL::nullPtr ) {
      krylovName_ = Glist.sublist("Krylov").get("Type","Conjugate Gradients");
      ekv_ = StringToEKrylov(krylovName_);
      krylov_ = KrylovFactory<Real>(parlist);
    }
    else {
      krylovName_ =  Glist.sublist("Krylov").get("User Defined Krylov Name",
                                                 "Unspecified User Defined Krylov Method");
    }
  }

  void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state ) {
    Step<Real>::initialize(x,s,g,obj,bnd,algo_state);
    if ( useSecantPrecond_ ) {
      gp_ = g.clone();
    }
  }

  void compute( Vector<Real> &s, const Vector<Real> &x,
                Objective<Real> &obj, BoundConstraint<Real> &bnd,
                AlgorithmState<Real> &algo_state ) {
    Real one(1);
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    // Build Hessian and Preconditioner object
    ROL::Ptr<Objective<Real> > obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<LinearOperator<Real> > hessian
      = ROL::makePtr<HessianNK>(obj_ptr,algo_state.iterateVec);
    ROL::Ptr<LinearOperator<Real> > precond;
    if ( useSecantPrecond_ ) {
      precond = secant_;
    }
    else {
      precond = ROL::makePtr<PrecondNK>(obj_ptr,algo_state.iterateVec);
    }

    // Run Krylov method
    flagKrylov_ = 0;
    krylov_->run(s,*hessian,*(step_state->gradientVec),*precond,iterKrylov_,flagKrylov_);

    // Check Krylov flags
    if ( flagKrylov_ == 2 && iterKrylov_ <= 1 ) {
      s.set((step_state->gradientVec)->dual());
    }
    s.scale(-one);
  }

  void update( Vector<Real> &x, const Vector<Real> &s,
               Objective<Real> &obj, BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    step_state->SPiter = iterKrylov_;
    step_state->SPflag = flagKrylov_;

    // Update iterate
    algo_state.iter++;
    x.plus(s);
    (step_state->descentVec)->set(s);
    algo_state.snorm = s.norm();

    // Compute new gradient
    if ( useSecantPrecond_ ) {
      gp_->set(*(step_state->gradientVec));
    }
    obj.update(x,true,algo_state.iter);
    if ( computeObj_ ) {
      algo_state.value = obj.value(x,tol);
      algo_state.nfval++;
    }
    obj.gradient(*(step_state->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update Secant Information
    if ( useSecantPrecond_ ) {
      secant_->updateStorage(x,*(step_state->gradientVec),*gp_,s,algo_state.snorm,algo_state.iter+1);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = step_state->gradientVec->norm();
  }

  std::string printHeader( void ) const {
    std::stringstream hist;

    if( verbosity_>0 ) {
      hist << std::string(109,'-') <<  "\n";
      hist << EDescentToString(DESCENT_NEWTONKRYLOV);
      hist << " status output definitions\n\n";
      hist << "  iter     - Number of iterates (steps taken) \n";
      hist << "  value    - Objective function value \n";
      hist << "  gnorm    - Norm of the gradient\n";
      hist << "  snorm    - Norm of the step (update to optimization vector)\n";
      hist << "  #fval    - Cumulative number of times the objective function was evaluated\n";
      hist << "  #grad    - Number of times the gradient was computed\n";
      hist << "  iterCG   - Number of Krylov iterations used to compute search direction\n";
      hist << "  flagCG   - Krylov solver flag" << "\n";
      hist << std::string(109,'-') << "\n";
    }

    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
    hist << "\n";
    return hist.str();
  }
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << EDescentToString(DESCENT_NEWTONKRYLOV);
    hist << " using " << krylovName_;
    if ( useSecantPrecond_ ) { 
      hist << " with " << ESecantToString(esec_) << " preconditioning";
    } 
    hist << "\n";
    return hist.str();
  }
  std::string print( AlgorithmState<Real> &algo_state, bool print_header = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( print_header ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::setw(10) << std::left << algo_state.nfval;
      hist << std::setw(10) << std::left << algo_state.ngrad;
      hist << std::setw(10) << std::left << iterKrylov_;
      hist << std::setw(10) << std::left << flagKrylov_;
      hist << "\n";
    }
    return hist.str();
  }
}; // class NewtonKrylovStep

} // namespace ROL

#endif
