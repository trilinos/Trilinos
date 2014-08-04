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

#ifndef ROL_TRUSTREGIONSTEP_H
#define ROL_TRUSTREGIONSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_TrustRegion.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::TrustRegionStep
    \brief Provides the interface to compute optimization steps
           with trust regions.
*/


namespace ROL {

template <class Real>
class TrustRegionStep : public Step<Real> {
private:

  Teuchos::RCP<Secant<Real> >       secant_;
  Teuchos::RCP<TrustRegion<Real> >  trustRegion_;

  ETrustRegion      etr_;        // Trust-Region Subproblem Solver Type
  ESecant           esec_;       // Secant Type

  bool useSecantHessVec_;
  bool useSecantPrecond_;

  bool useProjectedGrad_;

  std::vector<bool> useInexact_; // Inexactness Information
  int               TRflag_  ;   // Trust-Region Exit Flag
  int               TR_nfval_;   // Trust-Region Function Evaluation Number
  int               TR_ngrad_;   // Trust-Region Gradient Evaluation Number
  int               CGflag_;     // CG Termination Flag
  int               CGiter_;     // CG Iteration Count

  Real              alpha_init_; // Initial Line Search Parameter for Projected Methods
  int               max_fval_;   // Maximum Function Evaluations for Line Search              

  Real              scale0_;
  Real              scale1_;

  void updateGradient( Vector<Real> &x, Objective<Real> &obj, Constraints<Real> &con, 
                       AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    if ( this->useInexact_[1] ) {
      Real c = this->scale0_*std::max(1.e-2,std::min(1.0,1.e4*algo_state.gnorm));
      Real gtol1  = c*(state->searchSize);
      Real gtol0  = this->scale1_*gtol1 + 1.0;
      while ( gtol0 > gtol1*this->scale1_ ) {
        obj.gradient(*(state->gradientVec),x,gtol1);
        algo_state.gnorm = this->computeCriticalityMeasure(*(state->gradientVec),x,con);
        gtol0 = gtol1;
        c = this->scale0_*std::max(1.e-2,std::min(1.0,1.e4*algo_state.gnorm));
        gtol1 = c*std::min(algo_state.gnorm,state->searchSize);
      }
      algo_state.ngrad++;
    }
    else {
      Real gtol = std::sqrt(ROL_EPSILON);
      obj.gradient(*(state->gradientVec),x,gtol);
      algo_state.ngrad++;
      algo_state.gnorm = this->computeCriticalityMeasure(*(state->gradientVec),x,con);
    }
  }

  Real computeCriticalityMeasure( const Vector<Real> &g, const Vector<Real> &x, Constraints<Real> &con ) {
    if ( con.isActivated() ) {
      Teuchos::RCP<Vector<Real> > xnew = x.clone();
      if ( this->useProjectedGrad_ ) {
        xnew->set(g);
        con.computeProjectedGradient( *xnew, x );
      }
      else {
        xnew->set(x);
        xnew->axpy(-1.0,g);
        con.project(*xnew);
        xnew->axpy(-1.0,x);
      }
      return xnew->norm();
    }
    else {
      return g.norm();
    }
  }

public:

  virtual ~TrustRegionStep() {}

  TrustRegionStep( Teuchos::ParameterList & parlist ) : Step<Real>() {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Enumerations
    etr_   = StringToETrustRegion(parlist.get("Trust-Region Subproblem Solver Type","Cauchy Point"));  
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    // Secant Information
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    // Trust-Region Parameters
    step_state->searchSize = parlist.get("Initial Trust-Region Radius", -1.0);
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
    this->scale0_ = parlist.get("Gradient Update Tolerance Scaling",1.e-1);
    this->scale1_ = parlist.get("Gradient Update Relative Tolerance",2.0);     
     
    // Initialize Trust Region Subproblem Solver Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    max_fval_         = parlist.get("Maximum Number of Function Evaluations", 20);
    alpha_init_       = parlist.get("Initial Linesearch Parameter", 1.0);
    trustRegion_      = Teuchos::rcp( new TrustRegion<Real>(parlist) );

    // Secant Parameters
    secant_ = Teuchos::null;
    if ( useSecantPrecond_ || useSecantHessVec_ ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
  }

  TrustRegionStep( Teuchos::RCP<Secant<Real> > &secant, Teuchos::ParameterList &parlist ) 
    : Step<Real>(), secant_(secant) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Enumerations
    etr_   = StringToETrustRegion(parlist.get("Trust-Region Subproblem Solver Type","Cauchy Point"));  
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    // Secant Information
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    // Trust-Region Parameters
    step_state->searchSize = parlist.get("Initial Trust-Region Radius", -1.0);
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
    this->scale0_ = parlist.get("Gradient Update Tolerance Scaling",1.e-1);
    this->scale1_ = parlist.get("Gradient Update Relative Tolerance",2.0);     

    // Initialize Trust Region Subproblem Solver Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    max_fval_         = parlist.get("Maximum Number of Function Evaluations", 20);
    alpha_init_       = parlist.get("Initial Linesearch Parameter", 1.0);
    trustRegion_      = Teuchos::rcp( new TrustRegion<Real>(parlist) );
  }

  /** \brief Initialize step.
  */
  void initialize( Vector<Real> &x, Objective<Real> &obj, Constraints<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    algo_state.nfval = 0;
    algo_state.ngrad = 0;

    Real htol = std::sqrt(ROL_EPSILON);
    Real ftol = ROL_OVERFLOW; 

    step_state->descentVec  = x.clone();
    step_state->gradientVec = x.clone();

    if ( con.isActivated() ) {
      con.project(x);
    }

    // Update approximate gradient and approximate objective function.
    obj.update(x,true,algo_state.iter);    
    this->updateGradient(x,obj,con,algo_state);
    algo_state.snorm = 1.e10;
    algo_state.value = obj.value(x,ftol); 
    algo_state.nfval++;

    // Evaluate Objective Function at Cauchy Point
    if ( step_state->searchSize <= 0.0 ) {
      Teuchos::RCP<Vector<Real> > Bg = x.clone();
      if ( this->useSecantHessVec_ ) {
        this->secant_->applyB(*Bg,*(step_state->gradientVec),x);
      }
      else {
        obj.hessVec(*Bg,*(step_state->gradientVec),x,htol);
      }
      Real gBg = Bg->dot(*(step_state->gradientVec));
      Real alpha = 1.0;
      if ( gBg > ROL_EPSILON ) {
        alpha = algo_state.gnorm*algo_state.gnorm/gBg;
      }
      // Evaluate the objective function at the Cauchy point
      Teuchos::RCP<Vector<Real> > cp = x.clone();
      cp->set(*(step_state->gradientVec)); 
      cp->scale(-alpha);
      Teuchos::RCP<Vector<Real> > xnew = x.clone();
      xnew->set(x);
      xnew->plus(*cp);
      obj.update(*xnew);
      Real fnew = obj.value(*xnew,ftol); // MUST DO SOMETHING HERE WITH FTOL
      algo_state.nfval++;
      // Perform cubic interpolation to determine initial trust region radius
      Real gs = (step_state->gradientVec)->dot(*cp);
      Real a  = fnew - algo_state.value - gs - 0.5*alpha*alpha*gBg;
      if ( std::abs(a) < ROL_EPSILON ) { 
        // a = 0 implies the objective is quadratic in the negative gradient direction
        step_state->searchSize = alpha*algo_state.gnorm;
      }
      else {
        Real b  = 0.5*alpha*alpha*gBg;
        Real c  = gs;
        if ( b*b-3.0*a*c > ROL_EPSILON ) {
          // There is at least one critical point
          Real t1 = (-b-std::sqrt(b*b-3.0*a*c))/(3.0*a);
          Real t2 = (-b+std::sqrt(b*b-3.0*a*c))/(3.0*a);
          if ( 6.0*a*t1 + 2.0*b > 0.0 ) {
            // t1 is the minimizer
            step_state->searchSize = t1*alpha*algo_state.gnorm;          
          }
          else {
            // t2 is the minimizer
            step_state->searchSize = t2*alpha*algo_state.gnorm;
          }
        }
        else {
          step_state->searchSize = alpha*algo_state.gnorm;
        }
      }
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, Constraints<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    ProjectedObjective<Real> pObj(obj,con,this->secant_,this->useSecantPrecond_,this->useSecantHessVec_,eps);

    this->CGflag_ = 0;
    this->CGiter_ = 0;
    this->trustRegion_->run(s,algo_state.snorm,step_state->searchSize,this->CGflag_,this->CGiter_,
                            x,*(step_state->gradientVec),algo_state.gnorm,pObj);
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, Constraints<Real> &con, 
               AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();

    Real tol = std::sqrt(ROL_EPSILON);

    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    ProjectedObjective<Real> pObj(obj,con,this->secant_,this->useSecantPrecond_,this->useSecantHessVec_,eps);

    // Store previous step for constraint computations
    Teuchos::RCP<Vector<Real> > xold;
    if ( con.isActivated() ) {
      xold = x.clone();
      xold->set(x);
    }

    // Update trust-region information
    this->TRflag_   = 0;
    this->TR_nfval_ = 0;
    this->TR_ngrad_ = 0;
    Real fold = algo_state.value;
    Real fnew = 0.0;
    this->trustRegion_->update(x,fnew,state->searchSize,this->TR_nfval_,this->TR_ngrad_,this->TRflag_,
                               s,algo_state.snorm,fold,*(state->gradientVec),pObj);
    algo_state.value = fnew;
    algo_state.nfval += this->TR_nfval_;
    algo_state.ngrad += this->TR_ngrad_;
    algo_state.iter++;

    // Compute new gradient and update secant storage
    Teuchos::RCP<Vector<Real> > gp;
    if ( this->TRflag_ == 0 || this->TRflag_ == 1 ) {  
      // Perform line search (smoothing) to ensure decrease 
      if ( con.isActivated() ) {
        // Compute new gradient
        gp = x.clone(); 
        obj.gradient(*gp,x,tol); // MUST DO SOMETHING HERE WITH TOL
        algo_state.ngrad++;
        // Compute smoothed step
        Real alpha = 1.0;
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        xnew->set(x);
        xnew->axpy(-alpha*this->alpha_init_,*gp);
        con.project(*xnew);
        // Compute new objective value
        Real ftmp = obj.value(*xnew,tol); // MUST DO SOMETHING HERE WITH TOL
        algo_state.nfval++;
        // Perform smoothing
        int cnt = 0;
        alpha = 1.0/this->alpha_init_;
        while ( (fnew-ftmp) <= 1.e-4*(fnew-fold) ) { 
          xnew->set(x);
          xnew->axpy(-alpha*this->alpha_init_,*gp);
          con.project(*xnew);
          ftmp = obj.value(*xnew,tol); // MUST DO SOMETHING HERE WITH TOL
          algo_state.nfval++;
          if ( cnt >= this->max_fval_ ) {
            break;
          }
          alpha *= 0.5;
          cnt++;
        }
        // Store objective function and iteration information
        fnew = ftmp;
        x.set(*xnew);
      }

      // Store previous gradient for secant update
      if ( this->secant_ != Teuchos::null ) {
        gp = x.clone();
        gp->set(*(state->gradientVec));
      }

      // Update objective function and approximate model
      obj.update(x,true,algo_state.iter);
      this->updateGradient(x,obj,con,algo_state);

      // Update secant information
      if ( this->secant_ != Teuchos::null ) {
        if ( con.isActivated() ) { // Compute new constrained step
          Teuchos::RCP<Vector<Real> > st;
          st->set(x);
          st->axpy(-1.0,*xold);
          secant_->update(*(state->gradientVec),*gp,*st,algo_state.snorm,algo_state.iter+1);
        }
        else {
          secant_->update(*(state->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
        }
      }
    }    
  
    // Update algorithm state
    (algo_state.iterateVec)->set(x);
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(15) << std::left << "delta";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "tr_flag";
    if ( this->etr_ == TRUSTREGION_TRUNCATEDCG ) {
      hist << std::setw(10) << std::left << "iterCG";
      hist << std::setw(10) << std::left << "flagCG";
    }
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << ETrustRegionToString(this->etr_) << " Trust-Region solver";
    if ( this->useSecantPrecond_ || this->useSecantHessVec_ ) {
      if ( this->useSecantPrecond_ && !this->useSecantHessVec_ ) {
        hist << " with " << ESecantToString(this->esec_) << " preconditioning\n";
      }
      else if ( !this->useSecantPrecond_ && this->useSecantHessVec_ ) {
        hist << " with " << ESecantToString(this->esec_) << " Hessian approximation\n";
      }
      else {
        hist << " with " << ESecantToString(this->esec_) << " preconditioning and Hessian approximation\n";
      }
    }
    else {
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    const Teuchos::RCP<const StepState<Real> >& step_state = Step<Real>::getStepState();

    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << this->printName();
    }
    if ( printHeader ) {
      hist << this->printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << " "; 
      hist << std::setw(15) << std::left << step_state->searchSize; 
      hist << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(6)  << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(15) << std::left << step_state->searchSize; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << this->TRflag_;              
      if ( this->etr_ == TRUSTREGION_TRUNCATEDCG ) {
        hist << std::setw(10) << std::left << this->CGiter_;
        hist << std::setw(10) << std::left << this->CGflag_;
      }
      hist << "\n";
    }
    return hist.str();
  }

}; // class Step

} // namespace ROL

#endif
