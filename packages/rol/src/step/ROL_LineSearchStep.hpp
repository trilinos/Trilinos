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

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_Krylov.hpp"
#include "ROL_NonlinearCG.hpp"
#include "ROL_LineSearch.hpp"

#include <sstream>
#include <iomanip>

/** \class ROL::LineSearchStep
    \brief Provides the interface to compute optimization steps
           with line search.
*/


namespace ROL {

template <class Real>
class LineSearchStep : public Step<Real> {
private:

  Teuchos::RCP<Secant<Real> >              secant_;
  Teuchos::RCP<Krylov<Real> >              krylov_;
  Teuchos::RCP<NonlinearCG<Real> >         nlcg_;
  Teuchos::RCP<LineSearch<Real> >          lineSearch_;

  int iterKrylov_;
  int flagKrylov_;

  ELineSearch         els_;
  ENonlinearCG        enlcg_; 
  ECurvatureCondition econd_;
  EDescent            edesc_;
  ESecant             esec_;
  EKrylov             ekv_;
 
  int ls_nfval_;
  int ls_ngrad_;

  bool useSecantHessVec_;
  bool useSecantPrecond_;

  bool useProjectedGrad_;

  std::vector<bool> useInexact_;

public:

  virtual ~LineSearchStep() {}

  LineSearchStep( Teuchos::ParameterList &parlist ) : Step<Real>() {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    lineSearch_ = Teuchos::rcp( new LineSearch<Real>(parlist) );

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      Real CGtol1 = parlist.get("Absolute Krylov Tolerance", 1.e-4);
      Real CGtol2 = parlist.get("Relative Krylov Tolerance", 1.e-2);
      int maxitCG = parlist.get("Maximum Number of Krylov Iterations", 20);
      krylov_ = Teuchos::rcp( new Krylov<Real>(ekv_,CGtol1,CGtol2,maxitCG,useInexact_[2]) );
    }

    // Initialize Secant Object
    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }
  }

  LineSearchStep( Teuchos::RCP<Secant<Real> > &secant, Teuchos::ParameterList &parlist ) 
    : Step<Real>(), secant_(secant) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));

    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    lineSearch_ = Teuchos::rcp( new LineSearch<Real>(parlist) );

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioner", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      Real CGtol1 = parlist.get("Absolute Krylov Tolerance", 1.e-4);
      Real CGtol2 = parlist.get("Relative Krylov Tolerance", 1.e-2);
      int maxitCG = parlist.get("Maximum Number of Krylov Iterations", 20);
      krylov_ = Teuchos::rcp( new Krylov<Real>(ekv_,CGtol1,CGtol2,maxitCG,useInexact_[2]) );
    }

    // Secant Information
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }
     
    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    Real tol = std::sqrt(ROL_EPSILON);
    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    ProjectedObjective<Real> pObj(obj,con,this->secant_,this->useSecantPrecond_,this->useSecantHessVec_,eps);

    // Compute step s
    switch(this->edesc_) {
      case DESCENT_NEWTONKRYLOV:
        this->flagKrylov_ = 0;
        this->krylov_->run(s,this->iterKrylov_,this->flagKrylov_,*(step_state->gradientVec),x,pObj);
        break;
      case DESCENT_NEWTON:
      case DESCENT_SECANT:
        pObj.reducedInvHessVec(s,*(step_state->gradientVec),x,*(step_state->gradientVec),x,tol);
        break;
      case DESCENT_NONLINEARCG:
        this->nlcg_->run(s,*(step_state->gradientVec),x,obj);
        break;
      case DESCENT_STEEPEST:
        s.set(*(step_state->gradientVec));
        break;
      default: break;
    }

    // Compute g.dot(s)
    Real gs = 0.0;
    if ( !con.isActivated() ) {
      gs = -(step_state->gradientVec)->dot(s);
    }
    else {
      Teuchos::RCP<Vector<Real> > d = x.clone();
      if ( this->edesc_ == DESCENT_STEEPEST ) {
        d->set(x);
        d->axpy(-1.0,s);
        con.project(*d);
        d->scale(-1.0);
        d->plus(x);
        //d->set(s);
        //con.pruneActive(*d,s,x,eps);
        //con.pruneActive(*d,*(step_state->gradientVec),x,eps);
        gs = -(step_state->gradientVec)->dot(*d);
      }
      else {
        d->set(s);
        con.pruneActive(*d,*(step_state->gradientVec),x,eps);
        gs = -(step_state->gradientVec)->dot(*d);
        d->set(x);
        d->axpy(-1.0,*(step_state->gradientVec));
        con.project(*d);
        d->scale(-1.0);
        d->plus(x);
        con.pruneInactive(*d,*(step_state->gradientVec),x,eps);
        gs -= (step_state->gradientVec)->dot(*d);
      }
    }
    this->lineSearch_->setData((step_state->gradientVec),eps);

    // Check if s is a descent direction i.e., g.dot(s) < 0
    if ( gs >= 0.0 || (this->flagKrylov_ == 2 && this->iterKrylov_ <= 1) ) {
      s.set(*(step_state->gradientVec));
      if ( con.isActivated() ) {
        Teuchos::RCP<Vector<Real> > d = x.clone();
        d->set(s);
        con.pruneActive(*d,s,x);
        gs = -(step_state->gradientVec)->dot(*d);
      }
      else {
        gs = -(step_state->gradientVec)->dot(s);
      }
    }
    s.scale(-1.0);

    // Perform line search
    Real fnew  = algo_state.value;
    this->ls_nfval_ = 0;
    this->ls_ngrad_ = 0;
    this->lineSearch_->run(step_state->searchSize,fnew,this->ls_nfval_,this->ls_ngrad_,gs,s,x,obj,con);
    algo_state.nfval += this->ls_nfval_;
    algo_state.ngrad += this->ls_ngrad_;

    // Compute get scaled descent direction
    s.scale(step_state->searchSize);
    if ( con.isActivated() ) {
      s.plus(x);
      con.project(s);
      s.axpy(-1.0,x);
    }

    // Update step state information
    (step_state->descentVec)->set(s);

    // Update algorithm state information
    algo_state.snorm = s.norm();
    algo_state.value = fnew;
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON);
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Update iterate
    algo_state.iter++;
    x.axpy(1.0, s);
    obj.update(x,true,algo_state.iter);

    // Compute new gradient
    Teuchos::RCP<Vector<Real> > gp;
    if ( this->edesc_ == DESCENT_SECANT || (this->edesc_ == DESCENT_NEWTONKRYLOV && this->useSecantPrecond_) ) {
      gp = x.clone();
      gp->set(*(step_state->gradientVec));
    }
    obj.gradient(*(step_state->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update Secant Information
    if ( this->edesc_ == DESCENT_SECANT || (this->edesc_ == DESCENT_NEWTONKRYLOV && this->useSecantPrecond_) ) {
      secant_->update(*(step_state->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    if ( con.isActivated() ) {
      if ( this->useProjectedGrad_ ) {
        Teuchos::RCP<Vector<Real> > pg = x.clone();
        pg->set(*(step_state->gradientVec));
        con.computeProjectedGradient( *pg, x );
        algo_state.gnorm = pg->norm();
      }
      else {
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        xnew->set(x);
        xnew->axpy(-1.0,*(step_state->gradientVec));
        con.project(*xnew);
        xnew->axpy(-1.0,x);
        algo_state.gnorm = xnew->norm();
      }
    }
    else {
      algo_state.gnorm = (step_state->gradientVec)->norm();
    }
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6) << std::left << "iter";  
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm"; 
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "ls_#fval";
    hist << std::setw(10) << std::left << "ls_#grad";
    if ( this->edesc_ == DESCENT_NEWTONKRYLOV ) {
      hist << std::setw(10) << std::left << "iterCG";
      hist << std::setw(10) << std::left << "flagCG";
    }
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << EDescentToString(this->edesc_) 
         << " with " << ELineSearchToString(this->els_) 
         << " Linesearch satisfying " 
         << ECurvatureConditionToString(this->econd_) << "\n";
    if ( this->edesc_ == DESCENT_NEWTONKRYLOV ) {
      hist << "Krylov Type: " << EKrylovToString(this->ekv_) << "\n";
    }
    if ( this->edesc_ == DESCENT_SECANT || 
        (this->edesc_ == DESCENT_NEWTONKRYLOV && (this->useSecantPrecond_ || this->useSecantHessVec_)) ) {
      hist << "Secant Type: " << ESecantToString(this->esec_) << "\n";
    }
    if ( this->edesc_ == DESCENT_NONLINEARCG ) {
      hist << "Nonlinear CG Type: " << ENonlinearCGToString(this->enlcg_) << "\n";
    }
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
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
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(6) << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << this->ls_nfval_;              
      hist << std::setw(10) << std::left << this->ls_ngrad_;              
      if ( this->edesc_ == DESCENT_NEWTONKRYLOV ) {
        hist << std::setw(10) << std::left << this->iterKrylov_;
        hist << std::setw(10) << std::left << this->flagKrylov_;
      }
      hist << "\n";
    }
    return hist.str();
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
