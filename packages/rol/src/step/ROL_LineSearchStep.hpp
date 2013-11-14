//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_Krylov.hpp"
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

  Teuchos::RCP<Secant<Real> >     secant_;
  Teuchos::RCP<Krylov<Real> >     krylov_;
  Teuchos::RCP<LineSearch<Real> > lineSearch_;

  int iterKrylov_;
  int flagKrylov_;

  ELineSearch         els_; 
  ECurvatureCondition econd_;
  EDescent            edesc_;
  ESecant             esec_;
 
  int ls_nfval_;
  int ls_ngrad_;

public:

  virtual ~LineSearchStep() {}

  LineSearchStep( ELineSearch          els   = LINESEARCH_BACKTRACKING,
                  ECurvatureCondition  econd = CURVATURECONDITION_WOLFE,
                  EDescent             edesc = DESCENT_SECANT,              
                  int maxit = 20, Real c1 = 1.e-4, Real c2 = 0.9, Real LStol = 1.e-8, Real rho = 0.5,
                  ESecant esec = SECANT_LBFGS, int L = 10, int BBtype = 1,        // Secant Parameters
                  Real CGtol1 = 1.e-4, Real CGtol2 = 1.e-2, int maxitCG = 100 )
    : els_(els), econd_(econd), edesc_(edesc), esec_(esec) {

    lineSearch_ = Teuchos::rcp( new LineSearch<Real>( els_, econd_, edesc_, maxit, c1, c2, LStol, rho ) );

    krylov_ = Teuchos::null;
    if ( edesc_ == DESCENT_NEWTONKRYLOV || edesc_ == DESCENT_SECANTPRECOND ) {
      krylov_ = Teuchos::rcp( new Krylov<Real>(CGtol1,CGtol2,maxitCG) );
      iterKrylov_ = 0;
      flagKrylov_ = 0;
    }

    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) {
      if ( esec_ == SECANT_LBFGS ) {
        secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
      }
      else if ( esec_ == SECANT_LDFP ) {
        secant_ = Teuchos::rcp( new lDFP<Real>(L) );
      }
      else if ( esec_ == SECANT_LSR1 ) {
        secant_ = Teuchos::rcp( new lSR1<Real>(L) );
      }
      else if ( esec_ == SECANT_BARZILAIBORWEIN ) {
        secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
      }
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    // Compute step s
    if ( edesc_ == DESCENT_NEWTONKRYLOV || edesc_ == DESCENT_SECANTPRECOND ) {
      flagKrylov_ = 0;
      krylov_->CG(s,iterKrylov_,flagKrylov_,*(Step<Real>::state_->gradientVec),x,obj,secant_);
    }
    else if ( edesc_ == DESCENT_NEWTON ) {
      obj.invHessVec(s,*(Step<Real>::state_->gradientVec),x);
    }
    else if ( edesc_ == DESCENT_SECANT ) {
      secant_->applyH(s,*(Step<Real>::state_->gradientVec),x);
    }

    // Check if s is a descent direction
    Real gs = (Step<Real>::state_->gradientVec)->dot(s);
    if ( gs <= 0.0 || flagKrylov_ == 2 || edesc_ == DESCENT_STEEPEST ) {
      s.set(*(Step<Real>::state_->gradientVec));
      gs = (Step<Real>::state_->gradientVec)->dot(s);
    }
    s.scale(-1.0);
    gs *= -1.0;

    // Perform line search
    Real alpha = 0.0;
    Real fnew  = algo_state.value;
    ls_nfval_  = 0;
    ls_ngrad_  = 0;
    lineSearch_->run(alpha,fnew,ls_nfval_,ls_ngrad_,gs,s,x,obj);
    algo_state.nfval += ls_nfval_;
    algo_state.ngrad += ls_ngrad_;
    s.scale(alpha);

    // Update step state information
    (Step<Real>::state_->descentVec)->set(s);

    // Update algorithm state information
    algo_state.snorm = s.norm();
    algo_state.value = fnew;
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    // Update iterate
    x.axpy(1.0, s);

    // Compute new gradient
    Teuchos::RCP<Vector<Real> > gp;
    if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) {
      gp = x.clone();
      gp->set(*(Step<Real>::state_->gradientVec));
    }
    obj.gradient(*(Step<Real>::state_->gradientVec),x);
    algo_state.ngrad++;

    // Update Secant Information
    if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) {
      secant_->update(*(Step<Real>::state_->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = (Step<Real>::state_->gradientVec)->norm();
    algo_state.iter++;
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    std::stringstream hist;
    if ( algo_state.iter == 0 ) {
      hist << "\n" << EDescentToString(edesc_) 
           << " with " << ELineSearchToString(els_) 
           << " Linesearch satisfying " 
           << ECurvatureConditionToString(econd_) << "\n";
      if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) {
        hist << "Secant Type: " << ESecantToString(esec_) << "\n";
      }
      hist << "  ";
      hist << std::setw(6) << std::left << "iter";  
      hist << std::setw(15) << std::left << "value";
      hist << std::setw(15) << std::left << "gnorm"; 
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(10) << std::left << "#fval";
      hist << std::setw(10) << std::left << "#grad";
      hist << std::setw(10) << std::left << "ls_#fval";
      hist << std::setw(10) << std::left << "ls_#grad";
      if (    edesc_ == DESCENT_NEWTONKRYLOV 
           || edesc_ == DESCENT_SECANTPRECOND ) {
        hist << std::setw(10) << std::left << "iterCG";
        hist << std::setw(10) << std::left << "flagCG";
      }
      hist << "\n";
      hist << std::setfill('-') << std::setw(110) << "-" << "\n";
      hist << std::setfill(' ') << "  ";
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
      hist << std::setw(10) << std::left << ls_nfval_;              
      hist << std::setw(10) << std::left << ls_ngrad_;              
      if (    edesc_ == DESCENT_NEWTONKRYLOV 
           || edesc_ == DESCENT_SECANTPRECOND ) {
        hist << std::setw(10) << std::left << iterKrylov_;
        hist << std::setw(10) << std::left << flagKrylov_;
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
